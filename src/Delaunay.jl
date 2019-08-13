module Delaunay

using Random
using LightGraphs
using Plots
using LinearAlgebra
using DataStructures
using Logging

function is_closer(points, p, a::Int, b::Int)
    (points[p, :] .- points[a, :]) .^ 2 < (points[p, :] .- points[b, :]) .^ 2
end

function is_aligned(points, a::Int, b::Int, c::Int)
    v1 = points[b,:] .- points[a,:]
    v2 = points[c,:] .- points[a,:]
    v1[1] * v2[2] == v1[2] * v2[1]
end

"True if c is left of a -> b where a b c refere to points in points"
function is_left_of(points, a::Int, b::Int, c::Int)
    start_to_arrival = points[b,:] - points[a, :]
    norm = [-start_to_arrival[2] start_to_arrival[1]]
    (norm * (points[c,:]-points[b,:]))[1] >= 0
end

function is_strictly_left_of(points, a::Int, b::Int, c::Int)
    start_to_arrival = points[b,:] - points[a, :]
    norm = [-start_to_arrival[2] start_to_arrival[1]]
    (norm * (points[c,:]-points[b,:]))[1] > 0
end

"True if c is right of a -> b where a b c refere to points in points"
function is_right_of(points, a::Int, b::Int, c::Int)
    is_left_of(points, b, a, c)
end

function is_strictly_right_of(points, a::Int, b::Int, c::Int)
    is_strictly_left_of(points, b, a, c)
end

"Compute the intersection between a + λ * dir_a and b + γ * dir_b"
function line_intersection(a::Array{<:Real,1}, dir_a::Array{<:Real,1}, b::Array{<:Real,1}, dir_b::Array{<:Real,1})
    λ = [-dir_a dir_b]\(a .- b)
    a .+ λ[1] .* dir_a
end

line_intersection(points::Array{<:Real,1}, a::Int, b::Int, c::Int, d::Int) = line_intersection(points[a,:], points[b,:], points[c,:], points[d,:])

function circumcenter(points, a, b, c)
    xa = points[a, 1]
    ya = points[a, 2]
    xb = points[b, 1]
    yb = points[b, 2]
    xc = points[c, 1]
    yc = points[c, 2]
    M = [
        yb-ya  yb-yc
        xa-xb  xc-xb
    ]
    B = 1/2 * [xa-xc;ya-yc]
    λ = M\B
    1/2 * (points[a,:] + points[b,:]) + λ[1] * [ya-yb; xb-xa]
end

"True if p is in circumcircle a b c. (a,b,c ccw)"
function is_in_circumcircle(points, a::Int, b::Int, c::Int, p::Int)
    det(
        [
            points[a,1] points[a,2] (points[a,1]^2+points[a,2]^2) 1
            points[b,1] points[b,2] (points[b,1]^2+points[b,2]^2) 1
            points[c,1] points[c,2] (points[c,1]^2+points[c,2]^2) 1
            points[p,1] points[p,2] (points[p,1]^2+points[p,2]^2) 1
        ]
    ) > 0
end

mutable struct Triangle
    v1::Int64
    v2::Int64
    v3::Int64
    t1::Triangle
    t2::Triangle
    t3::Triangle
    function Triangle(a::Int, b::Int)
        t = new(a, b, 0)
        next = new(b, a, 0)
        t.t1 = t.t2 = t.t3 = next
        next.t1 = next.t2 = next.t3 = t
        t
    end
    function Triangle(v1::Int, v2::Int, v3::Int, t1::Triangle, t2::Triangle, t3::Triangle)
        new(v1, v2, v3, t1, t2, t3)
    end
end

function Triangle(points, a::Int, b::Int, c::Int)
    v1, v2, v3 = sort([a b c], dims=1)
    if is_left_of(points, v1, v2, v3)
        v2, v3 = v3, v2
    end
    t1 = Triangle(v2, v3) # We want to order the point of t1 such as t1.v1 -> t1.v2 is clock wise for t
    t2 = Triangle(v3, v1)
    t3 = Triangle(v1, v2)
    t = Triangle(v1, v2, v3, t1, t2, t3)

    # since ti.v1 -> ti.v2 is always clockwise, ti.t3 is always t
    t1.t3 = t2.t3 = t3.t3 = t
    t1.t1 = t3.t2 = t2
    t1.t2 = t2.t1 = t3
    t2.t2 = t3.t1 = t1
    t
end

Base.show(io::IO, t::Triangle) = print(io, "Triangle($(t.v1), $(t.v2), $(t.v3))")

circumcenter(points, triangle::Triangle) = circumcenter(points, triangle.v1, triangle.v2, triangle.v3)

function is_hull(triangle)
    triangle.v3 == 0
end

"Returns the next triangle on the hull, clock wise"
function next(triangle)
    triangle.t1
end

"Set the next triangle on the hull, clock wise"
function next!(prev, next)
    if !is_hull(prev)
        error("$prev is not on hull.")
    end
    if !is_hull(next)
        error("$next is not on hull.")
    end
    prev.t1 = next
    next.t2 = prev
end

"Returns the previous triangle on the hull, clock wise"
function prev(triangle)
    triangle.t2
end

"Set the previous triangle on the hull, clock wise"
function prev!(next, prev)
    next!(prev, next)
end

function replace_comp!(triangle, neighbour, comp)
    if triangle == neighbour.t1
        neighbour.t1 = comp
    elseif triangle == neighbour.t2
        neighbour.t2 = comp
    else
        neighbour.t3 = comp
    end
    comp.t3 = neighbour
end

replace_comp_t1!(triangle, comp) = replace_comp!(triangle, triangle.t1, comp)
replace_comp_t2!(triangle, comp) = replace_comp!(triangle, triangle.t2, comp)
replace_comp_t3!(triangle, comp) = replace_comp!(triangle, triangle.t3, comp)


"Delete a right set of triangle. `to_be_deleted` must be in the order of which
candidate choosing algorithm crosses the triangles.
"
function delete_right_triangles!(map, to_be_deleted, upper_triangle)
    next_right = next(upper_triangle)
    next_next_right = next(next_right)
    @debug "Starting deletion routine on right."
    for triangle in to_be_deleted
        pop!(map.triangles, triangle)
        @debug "Processing $triangle"
        if triangle == next_right # first iteration of the loop
            @debug "First iteration."
            continue
        elseif is_hull(triangle)
            @debug "Is hull, should be last iteration."
            continue
        elseif triangle.t2 == next_right
            @debug "triangle.t2 is current next_right."
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            comp_t3 = Triangle(triangle.v2, triangle.v1)
            @debug "Created $comp_t1 and $comp_t3."
            push!(map.triangles, comp_t1)
            push!(map.triangles, comp_t3)
            next!(upper_triangle, comp_t1)
            next!(comp_t1, comp_t3)
            next!(comp_t3, next_next_right)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t3!(triangle, comp_t3)
            next_right = comp_t1
            next_next_right = comp_t3
            rem_edge!(map.delaunay, triangle.v1, triangle.v3)
        elseif triangle.t3 == next_right
            @debug "triangle.t3 is current next_right."
            comp_t2 = Triangle(triangle.v1, triangle.v3)
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            @debug "Created $comp_t1 and $comp_t2."
            push!(map.triangles, comp_t1)
            push!(map.triangles, comp_t2)
            next!(upper_triangle, comp_t2)
            next!(comp_t2, comp_t1)
            next!(comp_t1, next_next_right)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t2!(triangle, comp_t2)
            next_right = comp_t2
            next_next_right = comp_t1
            rem_edge!(map.delaunay, triangle.v1, triangle.v2)
        else
            error("I am removing t1 on a right triangle.")
        end
    end
    # We might have deleted the next right triangle, especially if it was on hull
    push!(map.triangles, next(upper_triangle))
    @debug "Right deletion routine is over."
end

"Delete a left set of triangle. `to_be_deleted` must be in the order of which
candidate choosing algorithm crosses the triangles.
"
function delete_left_triangles!(map, to_be_deleted, upper_triangle)
    prev_left = prev(upper_triangle)
    prev_prev_left = prev(prev_left)
    @debug "Starting deletion routine on left."
    for triangle in to_be_deleted
        pop!(map.triangles, triangle)
        @debug "Processing $triangle"
        if triangle == prev_left # first iteration of the loop
            @debug "First iteration."
            continue
        elseif is_hull(triangle)
            @debug "Is hull, should be last iteration."
            continue
        elseif triangle.t1 == prev_left
            @debug "triangle.t1 is current prev_left."
            comp_t3 = Triangle(triangle.v2, triangle.v1)
            comp_t2 = Triangle(triangle.v1, triangle.v3)
            @debug "Created $comp_t2 and $comp_t3."
            push!(map.triangles, comp_t2)
            push!(map.triangles, comp_t3)
            replace_comp_t2!(triangle, comp_t2)
            replace_comp_t3!(triangle, comp_t3)
            prev!(upper_triangle, comp_t2)
            prev!(comp_t2, comp_t3)
            prev!(comp_t3, prev_prev_left)
            prev_left = comp_t2
            prev_prev_left = comp_t3
            rem_edge!(map.delaunay, triangle.v2, triangle.v3)
        elseif triangle.t2 == prev_left
            @debug "triangle.t2 is current prev_left."
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            comp_t3 = Triangle(triangle.v2, triangle.v1)
            @debug "Created $comp_t1 and $comp_t3."
            push!(map.triangles, comp_t1)
            push!(map.triangles, comp_t3)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t3!(triangle, comp_t3)
            prev!(upper_triangle, comp_t3)
            prev!(comp_t3, comp_t1)
            prev!(comp_t1, prev_prev_left)
            prev_left = comp_t3
            prev_prev_left = comp_t1
            rem_edge!(map.delaunay, triangle.v1, triangle.v3)
        else # triangle.t3 == prev_left
            @debug "triangle.t3 is current prev_left."
            comp_t2 = Triangle(triangle.v1, triangle.v3)
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            @debug "Created $comp_t1 and $comp_t2."
            push!(map.triangles, comp_t1)
            push!(map.triangles, comp_t2)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t2!(triangle, comp_t2)
            prev!(upper_triangle, comp_t1)
            prev!(comp_t1, comp_t2)
            prev!(comp_t2, prev_prev_left)
            prev_left = comp_t1
            prev_prev_left = comp_t2
            rem_edge!(map.delaunay, triangle.v1, triangle.v2)
        end
    end
    # We might have deleted the next right triangle, especially if it was on hull
    push!(map.triangles, prev(upper_triangle))
    @debug "Left deletion routine is over."
end

mutable struct Map{T<:Real}
    delaunay::SimpleGraph
    delaunay_points::Array{T,2}
    voronoi::SimpleGraph
    voronoi_points::Array{T,2}
    triangles::Set{Triangle}
    width::Real
    height::Real
end

function find_base_tangente(points, triangle_left, triangle_right)
    @debug "Starting find_base_tangente routine."
    prev_triangle_left = prev(triangle_left)
    prev_triangle_right = prev(triangle_right)
    there_is_a_point_below_left_next = is_strictly_right_of(
        points, triangle_left.v1, triangle_right.v1, triangle_left.v2
    )
    there_is_a_point_below_left_prev = is_strictly_right_of(
        points, triangle_left.v1, triangle_right.v1, prev_triangle_left.v1
    )
    there_is_a_point_below_right_next = is_strictly_right_of(
        points, triangle_left.v1, triangle_right.v1, triangle_right.v2
    )
    there_is_a_point_below_right_prev = is_strictly_right_of(
        points, triangle_left.v1, triangle_right.v1, prev_triangle_right.v1
    )
    while there_is_a_point_below_left_next || there_is_a_point_below_left_prev || there_is_a_point_below_right_next || there_is_a_point_below_right_prev
        @debug "Current tangente" triangle_left.v1 triangle_right.v1
        @debug "Working on" prev_triangle_left triangle_left prev_triangle_right triangle_right
        if there_is_a_point_below_left_next
            triangle_left, prev_triangle_left = next(triangle_left), triangle_left
        elseif there_is_a_point_below_left_prev
            triangle_left = prev_triangle_left
            prev_triangle_left = prev(triangle_left)

        end
        if there_is_a_point_below_right_next
            triangle_right, prev_triangle_right = next(triangle_right), triangle_right
        elseif there_is_a_point_below_right_prev
            triangle_right = prev_triangle_right
            prev_triangle_right = prev(triangle_right)
        end
        there_is_a_point_below_left_next = is_strictly_right_of(
            points, triangle_left.v1, triangle_right.v1, triangle_left.v2
        )
        there_is_a_point_below_left_prev = is_strictly_right_of(
            points, triangle_left.v1, triangle_right.v1, prev_triangle_left.v1
        )
        there_is_a_point_below_right_next = is_strictly_right_of(
            points, triangle_left.v1, triangle_right.v1, triangle_right.v2
        )
        there_is_a_point_below_right_prev = is_strictly_right_of(
            points, triangle_left.v1, triangle_right.v1, prev_triangle_right.v1
        )
    end

    # Corner case : aligned points
    right_next_aligned = is_aligned(points, triangle_left.v1, triangle_right.v1, triangle_right.v2)
    while right_next_aligned && is_closer(points, triangle_left.v1, triangle_right.v2, triangle_right.v1)
        prev_triangle_right, triangle_right = triangle_right, next(triangle_right)
        right_next_aligned = is_aligned(points, triangle_left.v1, triangle_right.v1, triangle_right.v2)
    end
    right_prev_aligned = is_aligned(points, triangle_left.v1, triangle_right.v1, prev_triangle_right.v1)
    while right_prev_aligned && is_closer(points, triangle_left.v1, prev_triangle_right.v1, triangle_right.v1)
        triangle_right = prev_triangle_right
        prev_triangle_right = prev(triangle_right)
        right_prev_aligned = is_aligned(points, triangle_left.v1, triangle_right.v1, prev_triangle_right.v1)
    end
    left_next_aligned = is_aligned(points, triangle_right.v1, triangle_left.v1, triangle_left.v2)
    while left_next_aligned && is_closer(points, triangle_right.v1, triangle_left.v2, triangle_left.v1)
        prev_triangle_left, triangle_left = triangle_left, next(triangle_left)
        left_next_aligned = is_aligned(points, triangle_right.v1, triangle_left.v1, triangle_left.v2)
    end
    left_prev_aligned = is_aligned(points, triangle_right.v1, triangle_left.v1, prev_triangle_left.v1)
    while left_prev_aligned && is_closer(points, triangle_right.v1, prev_triangle_left.v1, triangle_left.v1)
        triangle_left = prev_triangle_left
        prev_triangle_left = prev(triangle_left)
        left_prev_aligned = is_aligned(points, triangle_right.v1, triangle_left.v1, prev_triangle_left.v1)
    end
    @debug "End of find_base_tangente routine."
    triangle_left, triangle_right
end


function _in_circumcircle_left(points, triangle, base_left, base_right)
    if is_hull(triangle)
        false
    # We need to find out which vertex among v1 and v3 is base_left
    elseif triangle.v1 == base_left
        is_in_circumcircle(points, base_left, base_right, triangle.v3, triangle.v2)
    elseif triangle.v3 == base_left
        is_in_circumcircle(points, base_left, base_right, triangle.v2, triangle.v1)
    else # i.e. triangle.v2 == base_left
        is_in_circumcircle(points, base_left, base_right, triangle.v1, triangle.v3)
    end
end

function _in_circumcircle_right(points, triangle, base_left, base_right)
    if is_hull(triangle)
        false
    # We need to find out which vertex among v2 and v3 is base_right
    elseif triangle.v2 == base_right
        is_in_circumcircle(points, base_left, base_right, triangle.v3, triangle.v1)
    elseif triangle.v3 == base_right
        is_in_circumcircle(points, base_left, base_right, triangle.v1, triangle.v2)
    else # i.e triangle.v1 == base_right
        is_in_circumcircle(points, base_left, base_right, triangle.v2, triangle.v3)
    end
end

function _is_above_base_left(triangle, points, base_left, base_right)
    if is_hull(triangle) && triangle.v1 == base_left
        is_left_of(points, base_left, base_right, triangle.v2)
    elseif is_hull(triangle) || triangle.v3 == base_left
        is_left_of(points, base_left, base_right, triangle.v1)
    elseif triangle.v2 == base_left
        is_left_of(points, base_left, base_right, triangle.v3)
    else # i.e. v1 is base_left
        is_left_of(points, base_left, base_right, triangle.v2)
    end
end

function _is_above_base_right(triangle, points, base_left, base_right)
    if is_hull(triangle) && triangle.v2 == base_right
        is_left_of(points, base_left, base_right, triangle.v1)
    elseif is_hull(triangle) || triangle.v3 == base_right
        is_left_of(points, base_left, base_right, triangle.v2)
    elseif triangle.v1 == base_right
        is_left_of(points, base_left, base_right, triangle.v3)
    else # i.e. v2 is base_right
        is_left_of(points, base_left, base_right, triangle.v1)
    end
end

function _next_left(triangle, base_left)
    if is_hull(triangle) || triangle.v1 == base_left
        triangle.t3
    elseif triangle.v2 == base_left
        triangle.t1
    else # i.e. v3 is base_left
        triangle.t2
    end
end

function _next_right(triangle, base_right)
    if is_hull(triangle) || triangle.v2 == base_right
        triangle.t3
    elseif triangle.v3 == base_right
        triangle.t1
    else # i.e. v1 is base_right
        triangle.t2
    end
end

@enum Mode left=1 right=2
function find_candidate(mode::Mode, m, first_triangle, base_left, base_right)
    _in_circumcircle = if mode == left
        triangle -> _in_circumcircle_left(m.delaunay_points, triangle, base_left, base_right)
    else
        triangle -> _in_circumcircle_right(m.delaunay_points, triangle, base_left, base_right)
    end
    _is_above_base = if mode == left
        triangle -> _is_above_base_left(triangle, m.delaunay_points, base_left, base_right)
    else
        triangle -> _is_above_base_right(triangle, m.delaunay_points, base_left, base_right)
    end
    _next = if mode == left
        triangle -> _next_left(triangle, base_left)
    else
        triangle -> _next_right(triangle, base_right)
    end

    old_candidates = Deque{Triangle}()
    triangle = first_triangle
    found = false
    above_base = _is_above_base(triangle)
    while !(triangle in old_candidates) && !found && above_base
        next_triangle = _next(triangle)
        push!(old_candidates, triangle)
        if above_base && !_in_circumcircle(next_triangle)
            found=true
        else
            triangle = next_triangle
            above_base = _is_above_base(triangle)
        end
    end
    (found ? triangle : nothing), old_candidates
end

function _choose_candidate(points, base_left, base_right, left_candidate, right_candidate)
    b, c = base_left, base_right
    a = d = 0
    if is_hull(left_candidate) && left_candidate.v1 == b # Never happens if the hull is correctly build
        a = left_candidate.v2
    elseif is_hull(left_candidate) && left_candidate.v2 == b
        a = left_candidate.v1
    elseif left_candidate.v2 == b
        a = left_candidate.v3
    elseif left_candidate.v1 == b
        a = left_candidate.v2
    else # left_candidate.v3 == b
        a = left_candidate.v1
    end
    if is_hull(right_candidate) && right_candidate.v2 == c  # Never happens if the hull is correctly build
        d = right_candidate.v1
    elseif is_hull(right_candidate) && right_candidate.v1 == c
        d = right_candidate.v2
    elseif right_candidate.v1 == c
        d = right_candidate.v3
    elseif right_candidate.v2 == c
        d = right_candidate.v1
    else # right_candidate.v3 == c
        d = right_candidate.v2
    end
    @debug "Have to choose between $left_candidate, $right_candidate with test is_in_circumcircle(points, $a, $b, $c, $d)"
    if is_in_circumcircle(points, a, b, c, d) # if right candidate in circumcircle(left_candidate, base_left, base_right)
        right_candidate
    else
        left_candidate
    end
end

function triangulate!(m, first=1, last=-1)
    if last == -1
        last = div(length(m.delaunay_points), 2)
    end
    if first == last # never happens
        error("Trying to triangulate a single point $first.")
    elseif last - first == 1
        @debug "Triangulate on 2 vertices : $first - $last."
        add_edge!(m.delaunay, first, last)
        t = Triangle(first, last)
        push!(m.triangles, t)
        push!(m.triangles, t.t3)
        @debug "Created $t."
        @debug to_geogebra(m)
        t.t3, t
    elseif last - first == 2
        @debug "Triangulate on 3 vertices : $first to $last."
        if !is_aligned(m.delaunay_points, first, first+1, last)
            add_edge!(m.delaunay, first, first+1)
            add_edge!(m.delaunay, first+1, last)
            add_edge!(m.delaunay, last, first)
            t = Triangle(m.delaunay_points, first, first+1, last)
            push!(m.triangles, t)
            push!(m.triangles, t.t1)
            push!(m.triangles, t.t2)
            push!(m.triangles, t.t3)
            @debug "Created $t."
            @debug to_geogebra(m)
            if t.v2 < t.v3
                t.t3, t.t2
            else
                t.t3, t.t1
            end
        else
            add_edge!(m.delaunay, first, first+1)
            add_edge!(m.delaunay, first+1, last)
            t_first = Triangle(first, first+1)
            t_last = Triangle(first+1, last)
            next!(t_first, t_last)
            next!(t_last, t_last.t3)
            next!(t_last.t3, t_first.t3)
            next!(t_first.t3, t_first)
            push!(m.triangles, t_first)
            push!(m.triangles, t_first.t3)
            push!(m.triangles, t_last)
            push!(m.triangles, t_last.t3)
            @debug "Created $t_first and $t_last."
            @debug to_geogebra(m)
            if m.delaunay_points[first, 1] < m.delaunay_points[first+1, 1] # vertical line
                t_first, t_last.t3
            else # horizontal line
                t_first.t3, t_last
            end
        end

    else
        mid = div(first + last, 2)
        triangle_left_left, triangle_left_right = triangulate!(m, first, mid)
        triangle_right_left, triangle_right_right = triangulate!(m, mid+1, last)

        @debug "Triangulate on $(last - first + 1) vertices : $first to $last."

        base_triangle_left, base_triangle_right = find_base_tangente(m.delaunay_points, triangle_left_right, triangle_right_left)
        base_left = base_triangle_left.v1
        base_right = base_triangle_right.v1

        triangle_left = prev(base_triangle_left)
        triangle_right = base_triangle_right

        # Reminder : for 2edges triangles, we force v1 -> v2 to be clockwise
        base_triangle = Triangle(base_right, base_left)
        push!(m.triangles, base_triangle)
        next!(base_triangle, base_triangle_left)
        prev!(base_triangle, prev(base_triangle_right))
        lower_triangle = base_triangle

        upper_triangle = base_triangle.t3
        push!(m.triangles, upper_triangle)
        next!(upper_triangle, triangle_right)
        prev!(upper_triangle, triangle_left)

        right_candidate = next(upper_triangle)
        left_candidate = prev(upper_triangle)

        base_triangle_edge = :t3

        add_edge!(m.delaunay, base_left, base_right)

        while !isnothing(right_candidate) || !isnothing(left_candidate)
            @debug "Base triangle is $base_triangle."
            right_candidate, delete_triangle_right = find_candidate(right, m, triangle_right, base_left, base_right)
            left_candidate, delete_triangle_left = find_candidate(left, m, triangle_left, base_left, base_right)

            @debug "Left candidate is " left_candidate
            @debug "Right candidate is " right_candidate

            chosen = :none

            to_be_deleted = isnothing(right_candidate) ? delete_triangle_left : delete_triangle_right
            triangle = isnothing(right_candidate) ? left_candidate : right_candidate

            if !isnothing(left_candidate) && !isnothing(right_candidate)
                triangle = _choose_candidate(m.delaunay_points, base_left, base_right, left_candidate, right_candidate)
                if triangle == left_candidate
                    to_be_deleted = delete_triangle_left
                end
            end
            if isnothing(triangle)
                chosen = :none
            elseif triangle == left_candidate
                chosen = :left
            elseif triangle == right_candidate
                chosen = :right
            end

            if chosen != :none
                @debug "Chose $chosen."

                new_base_left, new_base_right = base_left, base_right
                new_base_triangle_edge = :t1
                new_vertex = 0

                if chosen == :left
                    delete_left_triangles!(m, to_be_deleted, upper_triangle)
                    triangle = prev(upper_triangle)
                    # Now we are sure that `triangle` is on hull.
                    new_vertex = triangle.v1
                else
                    delete_right_triangles!(m, to_be_deleted, upper_triangle)
                    triangle = next(upper_triangle)
                    # Now we are sure that `triangle` is on hull.
                    new_vertex = triangle.v2
                end

                # Let's order the vertices of this new triangle
                triangle.v1, triangle.v2, triangle.v3 = sort([new_vertex, base_left, base_right])
                if is_left_of(m.delaunay_points, triangle.v1, triangle.v2, triangle.v3)
                    triangle.v2, triangle.v3 = triangle.v3, triangle.v2
                end

                # time to ensure the neighbors triangle are the right ones.

                upper_triangle.t3 = triangle
                if chosen == :left
                    # `prev(triangle)` is still valid since we have only
                    # modified its vertices for the moment.
                    new_left = prev(triangle)
                    if base_left == triangle.v1
                        upper_triangle.v1, upper_triangle.v2 = triangle.v2, triangle.v3
                        triangle.t1, triangle.t2 = upper_triangle, base_triangle
                        new_base_triangle_edge = :t1
                        new_base_left = triangle.v2
                    elseif base_left == triangle.v3
                        upper_triangle.v1, upper_triangle.v2 = triangle.v1, triangle.v2
                        triangle.t1, triangle.t2, triangle.t3 = base_triangle, triangle.t3, upper_triangle
                        new_base_triangle_edge = :t3
                        new_base_left = triangle.v1
                    end
                    prev!(upper_triangle, new_left)
                else
                    upper_triangle.v1, upper_triangle.v2 = triangle.v1, triangle.v2
                    next!(upper_triangle, next(triangle))
                    triangle.t1, triangle.t2, triangle.t3 = triangle.t3, base_triangle, upper_triangle
                    new_base_triangle_edge = :t3
                    new_base_right = triangle.v2
                end

                # Link base_triangle
                if base_triangle_edge == :t1
                    base_triangle.t1 = triangle
                elseif base_triangle_edge == :t2
                    base_triangle.t2 = triangle
                else
                    base_triangle.t3 = triangle
                end

                # We can finally add the new edge :D
                add_edge!(m.delaunay, new_base_left, new_base_right)

                # set up the new value for loop variables
                triangle_left = prev(upper_triangle)
                triangle_right = next(upper_triangle)
                base_triangle = triangle
                base_left, base_right = new_base_left, new_base_right
                base_triangle_edge = new_base_triangle_edge
            end
        end

        @debug to_geogebra(m)
        lower_triangle, upper_triangle
    end

end

function project_bissection_on_nearest_border(points::Array{<:Real,2}, triangle::Triangle, width::Real, height::Real)
    center = (points[triangle.v1,:] .+ points[triangle.v2,:]) ./ 2
    direction = [0 -1;1 0] * (points[triangle.v1,:] .- points[triangle.v2,:])
    sort([
            line_intersection(center, direction, [0;0], [0;height]),
            line_intersection(center, direction, [width;0], [0;height]),
            line_intersection(center, direction, [0;0], [width;0]),
            line_intersection(center, direction, [0;height], [width;0])
        ],
        lt=(x,y)-> norm(center .- x) < norm(center .- y)
    )[1]
end


function create_voronoi_from_triangles!(map::Map)
    processed = Set{Triangle}()
    stack = Stack{Triangle}()
    indices = Dict()
    for (i,triangle) in enumerate(map.triangles)
        push!(indices, triangle => i)
    end
    push!(stack, iterate(indices)[1].first)
    map.voronoi_points = zeros(Float64, length(indices), 2)
    map.voronoi = SimpleGraph(length(indices))
    while !isempty(stack)
        triangle = pop!(stack)
        push!(processed, triangle)
        index = get(indices, triangle, 0)
        map.voronoi_points[index,:] = if is_hull(triangle)
            project_bissection_on_nearest_border(map.delaunay_points, triangle, map.width, map.height)
        else
            circumcenter(map.delaunay_points, triangle)
        end

        if triangle.t1 ∉ processed
            if !(is_hull(triangle) && is_hull(triangle.t1))
                add_edge!(map.voronoi, index, get(indices, triangle.t1, 0))
            end
            push!(stack, triangle.t1)
        end
        if triangle.t2 ∉ processed
            if !(is_hull(triangle) && is_hull(triangle.t2))
                add_edge!(map.voronoi, index, get(indices, triangle.t2, 0))
            end
            push!(stack, triangle.t2)
        end
        if triangle.t3 ∉ processed
            if !(is_hull(triangle) && is_hull(triangle.t3))
                add_edge!(map.voronoi, index, get(indices, triangle.t3, 0))
            end
            push!(stack, triangle.t3)
        end
    end
end

function Map(points::Array{T,2}) where T <: Real
    m = Map(SimpleGraph(div(length(points),2)), sortslices(points, dims=1), SimpleGraph(), zeros(0,0), Set{Triangle}(), maximum(points[:,1]), maximum(points[:,2]))
    left,right = triangulate!(m)
    create_voronoi_from_triangles!(m)
    m
end

function Map(height::Number, width::Number, nb_points::Int; res::Float64=0.5)
    points = sortslices([rand(nb_points) .* width rand(nb_points) .* height], dims=1)
    @debug "Generating from" points
    delaunay = SimpleGraph(nb_points)
    m = Map(delaunay, points, SimpleGraph(), zeros(0,0), Set{Triangle}(), width, height)
    triangulate!(m)
    create_voronoi_from_triangles!(m)
    m
end


function plot_points(points, type=:delaunay; show_nums=true)
    plot()
    plot_points!(points, type, show_nums=show_nums)
end

function plot_points!(points, type=:delaunay; show_nums=true)
    marker = if type == :delaunay
        if show_nums
            (:circle, 20, 0.2, :orange)
        else
            (:circle, 100, 1, :red)
        end
    else
        if show_nums
            (:circle, 10, 0.2, :blue)
        else
            (:circle, 100, 1, :blue)
        end
    end
    if show_nums
        plot!(
            points[:,1], points[:,2],
            marker=marker, line=:scatter,
            series_annotations=collect(map(string, 1:length(points[:,1]))),
            aspect_ratio=:equal
        )
    else
        plot!(
            points[:,1], points[:,2],
            marker=marker, line=:scatter, aspect_ratio=:equal
        )
    end
end

function plot_map(m; show_nums=true, window=:fit)
    plot()
    plot_map!(m, show_nums=show_nums, window=window)
end
function plot_map!(m; show_nums=true, window=:fit)
    if window == :full
        plot!(
            aspect_ratio=:equal
        )
    else
        plot!(
            aspect_ratio=:equal,
            xlims=(0,m.width),
            ylims=(0,m.height),
        )
    end
    plot!(
        hcat(map(x->[m.delaunay_points[x.src, 1]; m.delaunay_points[x.dst, 1]],edges(m.delaunay))...),
        hcat(map(x->[m.delaunay_points[x.src, 2]; m.delaunay_points[x.dst, 2]],edges(m.delaunay))...),
        color="red",
        legend=false,
    )
    plot!(
        hcat(map(x->[m.voronoi_points[x.src, 1]; m.voronoi_points[x.dst, 1]],edges(m.voronoi))...),
        hcat(map(x->[m.voronoi_points[x.src, 2]; m.voronoi_points[x.dst, 2]],edges(m.voronoi))...),
        color="blue",
        legend=false,
    )
    plot_points!(m.voronoi_points, :voronoi, show_nums=show_nums)
    plot_points!(m.delaunay_points, show_nums=show_nums)
end

function to_geogebra(points)
    (*)(
        "Point[{",
        join(["{$x, $y}" for (x,y) in zip(points[:,1], points[:,2])], ", "),
        "}]"
    )
end

function to_geogebra(map::Map)
    convert_to_ascii(n) = Char(n + 64)
    n_delaunay = div(length(map.delaunay_points), 2)
    join(
        [
            to_geogebra(map.delaunay_points),
            (*)(
                "{",
                join(["Segment[$(convert_to_ascii(e.src)), $(convert_to_ascii(e.dst))]" for e in edges(map.delaunay)], ", "),
                "}"
            ),
            (*)(
                "{",
                join(["Polygone[$(convert_to_ascii(t.v1)), $(convert_to_ascii(t.v2)), $(convert_to_ascii(t.v3))]" for t in map.triangles if !is_hull(t)], ", "),
                "}"
            ),
            to_geogebra(map.voronoi_points),
            (*)(
                "{",
                join(["Segment[$(convert_to_ascii(e.src+n_delaunay)), $(convert_to_ascii(e.dst+n_delaunay))]" for e in edges(map.voronoi)], ", "),
                "}"
            ),
        ],
        "\n"
    )
end

"Dummy function to test if a Map is delaunay."
function is_delaunay(map)
    result = true
    for t in map.triangles
        if is_hull(t)
            continue
        end
        for p in 1:div(length(map.delaunay_points), 2)
            if !(p in (t.v1, t.v2, t.v3))
                ok = !is_in_circumcircle(map.delaunay_points, t.v1, t.v3, t.v2, p)
                if !ok
                    @debug "failed on triangle $t point $p"
                    result = false
                end
            end
        end
    end
    result
end

export Map, to_geogebra, plot_points, plot_points!, plot_map, plot_map!

end
