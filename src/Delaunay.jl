module Delaunay

using Random
using LightGraphs
using Plots
using LinearAlgebra
using DataStructures

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
function delete_right_triangles!(graph, to_be_deleted, upper_triangle)
    next_right = next(upper_triangle)
    next_next_right = next(next_right)
    for triangle in to_be_deleted
        if triangle == next_right # first iteration of the loop
            continue
        elseif is_hull(triangle)
            continue
        elseif triangle.t2 == next_right
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            comp_t3 = Triangle(triangle.v2, triangle.v1)
            next!(upper_triangle, comp_t1)
            next!(comp_t1, comp_t3)
            next!(comp_t3, next_next_right)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t3!(triangle, comp_t3)
            next_right = comp_t1
            next_next_right = comp_t3
            rem_edge!(graph, triangle.v1, triangle.v3)
        elseif triangle.t3 == next_right
            comp_t2 = Triangle(triangle.v1, triangle.v3)
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            next!(upper_triangle, comp_t2)
            next!(comp_t2, comp_t1)
            next!(comp_t1, next_next_right)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t2!(triangle, comp_t2)
            next_right = comp_t2
            next_next_right = comp_t1
            rem_edge!(graph, triangle.v1, triangle.v2)
        else
            error("I am removing t1 on a right triangle.")
        end
    end
end

"Delete a left set of triangle. `to_be_deleted` must be in the order of which
candidate choosing algorithm crosses the triangles.
"
function delete_left_triangles!(graph, to_be_deleted, upper_triangle)
    prev_left = prev(upper_triangle)
    prev_prev_left = prev(prev_left)
    for triangle in to_be_deleted
        if triangle == prev_left # first iteration of the loop
            continue
        elseif is_hull(triangle)
            continue
        elseif triangle.t1 == prev_left
            comp_t3 = Triangle(triangle.v2, triangle.v1)
            comp_t2 = Triangle(triangle.v1, triangle.v3)
            replace_comp_t2!(triangle, comp_t2)
            replace_comp_t3!(triangle, comp_t3)
            prev!(upper_triangle, comp_t2)
            prev!(comp_t2, comp_t3)
            prev!(comp_t3, prev_prev_left)
            prev_left = comp_t2
            prev_prev_left = comp_t3
            rem_edge!(graph, triangle.v2, triangle.v3)
        elseif triangle.t2 == prev_left
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            comp_t3 = Triangle(triangle.v2, triangle.v1)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t3!(triangle, comp_t3)
            prev!(upper_triangle, comp_t3)
            prev!(comp_t3, comp_t1)
            prev!(comp_t1, prev_prev_left)
            prev_left = comp_t3
            prev_prev_left = comp_t1
            rem_edge!(graph, triangle.v1, triangle.v3)
        else # triangle.t3 == next_right
            comp_t2 = Triangle(triangle.v1, triangle.v3)
            comp_t1 = Triangle(triangle.v3, triangle.v2)
            replace_comp_t1!(triangle, comp_t1)
            replace_comp_t2!(triangle, comp_t2)
            prev!(upper_triangle, comp_t1)
            prev!(comp_t1, comp_t2)
            prev!(comp_t2, prev_prev_left)
            prev_left = comp_t1
            prev_prev_left = comp_t2
            rem_edge!(graph, triangle.v1, triangle.v2)
        end
    end
end

struct Map{T<:Real}
    delaunay::SimpleGraph
    points::Array{T,2}
end

function find_base_tangente(points, triangle_left, triangle_right)
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
        triangle -> _in_circumcircle_left(m.points, triangle, base_left, base_right)
    else
        triangle -> _in_circumcircle_right(m.points, triangle, base_left, base_right)
    end
    _is_above_base = if mode == left
        triangle -> _is_above_base_left(triangle, m.points, base_left, base_right)
    else
        triangle -> _is_above_base_right(triangle, m.points, base_left, base_right)
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
    if is_hull(left_candidate) && left_candidate.v1 == b
        a = left_candidate.v2
    elseif left_candidate.v2 == b # also matches the cas of triangle on hull
        a = left_candidate.v1
    elseif left_candidate.v1 == b
        a = left_candidate.v3
    else # left_candidate.v3 == b
        a = left_candidate.v2
    end
    if is_hull(right_candidate) && right_candidate.v2 == c
        d = right_candidate.v1
    elseif right_candidate.v1 == c # also matches the cas of triangle on hull
        d = right_candidate.v2
    elseif right_candidate.v2 == c
        d = right_candidate.v3
    else # right_candidate.v3 == c
        d = right_candidate.v1
    end
    if is_in_circumcircle(points, a, b, c, d) # if right candidate in circumcircle(left_candidate, base_left, base_right)
        right_candidate
    else
        left_candidate
    end
end

function triangulate!(m, first=1, last=-1)
    if last == -1
        last = div(length(m.points), 2)
    end
    if first == last # never happens
        error("Trying to triangulate a single point $first.")
    elseif last - first == 1
        add_edge!(m.delaunay, first, last)
        t = Triangle(first, last)
        t.t3, t
    elseif last - first == 2
        if !is_aligned(m.points, first, first+1, last)
            add_edge!(m.delaunay, first, first+1)
            add_edge!(m.delaunay, first+1, last)
            add_edge!(m.delaunay, last, first)
            t = Triangle(m.points, first, first+1, last)
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
            if m.points[first, 1] < m.points[first+1, 1] # vertical line
                t_first, t_last.t3
            else # horizontal line
                t_first.t3, t_last
            end
        end

    else
        mid = div(first + last, 2)
        triangle_left_left, triangle_left_right = triangulate!(m, first, mid)
        triangle_right_left, triangle_right_right = triangulate!(m, mid+1, last)

        base_triangle_left, base_triangle_right = find_base_tangente(m.points, triangle_left_right, triangle_right_left)
        base_left = base_triangle_left.v1
        base_right = base_triangle_right.v1

        triangle_left = prev(base_triangle_left)
        triangle_right = base_triangle_right

        # Reminder : for 2edges triangles, we force v1 -> v2 to be clockwise
        base_triangle = Triangle(base_right, base_left)
        next!(base_triangle, base_triangle_left)
        prev!(base_triangle, prev(base_triangle_right))

        upper_triangle = base_triangle.t3
        next!(upper_triangle, triangle_right)
        prev!(upper_triangle, triangle_left)

        right_candidate = next(upper_triangle)
        left_candidate = prev(upper_triangle)

        base_triangle_edge = :t3

        add_edge!(m.delaunay, base_left, base_right)

        while !isnothing(right_candidate) || !isnothing(left_candidate)
            right_candidate, delete_triangle_right = find_candidate(right, m, triangle_right, base_left, base_right)
            left_candidate, delete_triangle_left = find_candidate(left, m, triangle_left, base_left, base_right)

            chosen = :none

            to_be_deleted = isnothing(right_candidate) ? delete_triangle_left : delete_triangle_right
            triangle = isnothing(right_candidate) ? left_candidate : right_candidate

            if !isnothing(left_candidate) && !isnothing(right_candidate)
                triangle = _choose_candidate(m.points, base_left, base_right, left_candidate, right_candidate)
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

                new_base_left, new_base_right = base_left, base_right
                new_base_triangle_edge = :t1
                new_vertex = 0

                if chosen == :left
                    delete_left_triangles!(m.delaunay, to_be_deleted, upper_triangle)
                    triangle = prev(upper_triangle)
                    # Now we are sure that `triangle` is on hull.
                    new_vertex = triangle.v1
                else
                    delete_right_triangles!(m.delaunay, to_be_deleted, upper_triangle)
                    triangle = next(upper_triangle)
                    # Now we are sure that `triangle` is on hull.
                    new_vertex = triangle.v2
                end

                # Let's order the vertices of this new triangle
                triangle.v1, triangle.v2, triangle.v3 = sort([new_vertex, base_left, base_right])
                if is_left_of(m.points, triangle.v1, triangle.v2, triangle.v3)
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


        while !is_hull(triangle_left_left)
            if is_hull(triangle_left_left.t1)
                triangle_left_left = triangle_left_left.t1
            elseif is_hull(triangle_left_left.t2)
                    triangle_left_left = triangle_left_left.t2
            elseif is_hull(triangle_left_left.t3)
                triangle_left_left = triangle_left_left.t3
            else
                triangle_left_left = triangle_left_left.t1
            end
        end
        while !is_hull(triangle_right_right)
            if is_hull(triangle_right_right.t1)
                triangle_right_right = triangle_right_right.t1
            elseif is_hull(triangle_right_right.t2)
                triangle_right_right = triangle_right_right.t2
            elseif is_hull(triangle_right_right.t3)
                triangle_right_right = triangle_right_right.t3
            else
                triangle_right_right = triangle_right_right.t2
            end
        end
        triangle_left_left, triangle_right_right
    end

end

function Map(points::Array{T,2}) where T <: Real
    m = Map(SimpleGraph(div(length(points),2)), sortslices(points, dims=1))
    left,right = triangulate!(m)
    m
end

function Map(height::Number, width::Number, nb_points::Int; res::Float64=0.5)
    points = sortslices([rand(nb_points) .* width rand(nb_points) .* height], dims=1)
    delaunay = SimpleGraph(nb_points)
    m = Map(delaunay, points)
    triangulate!(m)
    m
end


function plot_points(points)
    plot(
        points[:,1], points[:,2],
        marker=(20, 0.2, :orange), line=:scatter,
        series_annotations=collect(map(string, 1:length(points[:,1])))
    )
end

function plot_points!(points)
    plot!(
        points[:,1], points[:,2],
        marker=(20, 0.2, :orange), line=:scatter,
        series_annotations=collect(map(string, 1:length(points[:,1])))
    )
end

function plot_map(m)
    plot(
        hcat(map(x->[m.points[x.src, 1]; m.points[x.dst, 1]],edges(m.delaunay))...),
        hcat(map(x->[m.points[x.src, 2]; m.points[x.dst, 2]],edges(m.delaunay))...),
        color="red",
        legend=false
    )
    plot_points!(m.points)
end
function plot_map!(m)
    plot!(
        hcat(map(x->[m.points[x.src, 1]; m.points[x.dst, 1]],edges(m.delaunay))...),
        hcat(map(x->[m.points[x.src, 2]; m.points[x.dst, 2]],edges(m.delaunay))...),
        color="red",
        legend=false
    )
    plot_points!(m.points)
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
    join(
        [
            to_geogebra(map.points),
            (*)(
                "{",
                join(["Segment[$(convert_to_ascii(e.src)), $(convert_to_ascii(e.dst))]" for e in edges(map.delaunay)], ", "),
                "}"
            )
        ],
        "\n"
    )
end

export Map, to_geogebra

end
