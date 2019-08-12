using Test
using Delaunay

@testset "Tests for triangulation" begin
@testset "Test function is_closer" begin
    points = [
        0  0
        1  0
        1  1
    ]
    @test Delaunay.is_closer(points, 1, 2, 3) == true
    @test Delaunay.is_closer(points, 1, 3, 2) == false
end
@testset "Test function is_aligned" begin
    points = [
        0  0
        1  0
        1  1
        2  0
    ]
    @test Delaunay.is_aligned(points, 1, 2, 3) == false
    @test Delaunay.is_aligned(points, 1, 2, 4) == true
end
@testset "Test function is_left_of" begin
    points = [
         0  0
         1  0
         1  1
         1 -1
        -1  1
        -1 -1
         2  0
    ]
    @test Delaunay.is_left_of(points, 1, 2, 3) == true
    @test Delaunay.is_left_of(points, 1, 2, 4) == false
    @test Delaunay.is_left_of(points, 1, 2, 5) == true
    @test Delaunay.is_left_of(points, 1, 2, 6) == false
    @test Delaunay.is_left_of(points, 1, 2, 7) == true
end
@testset "Test function is_strictly_left_of" begin
    points = [
         0  0
         1  0
         1  1
         1 -1
        -1  1
        -1 -1
         2  0
    ]
    @test Delaunay.is_strictly_left_of(points, 1, 2, 3) == true
    @test Delaunay.is_strictly_left_of(points, 1, 2, 4) == false
    @test Delaunay.is_strictly_left_of(points, 1, 2, 5) == true
    @test Delaunay.is_strictly_left_of(points, 1, 2, 6) == false
    @test Delaunay.is_strictly_left_of(points, 1, 2, 7) == false
end
@testset "Test function is_right_of" begin
    points = [
         0  0
         1  0
         1  1
         1 -1
        -1  1
        -1 -1
         2  0
    ]
    @test Delaunay.is_right_of(points, 1, 2, 3) == false
    @test Delaunay.is_right_of(points, 1, 2, 4) == true
    @test Delaunay.is_right_of(points, 1, 2, 5) == false
    @test Delaunay.is_right_of(points, 1, 2, 6) == true
    @test Delaunay.is_right_of(points, 1, 2, 7) == true
end
@testset "Test function is_strictly_right_of" begin
    points = [
         0  0
         1  0
         1  1
         1 -1
        -1  1
        -1 -1
         2  0
    ]
    @test Delaunay.is_strictly_right_of(points, 1, 2, 3) == false
    @test Delaunay.is_strictly_right_of(points, 1, 2, 4) == true
    @test Delaunay.is_strictly_right_of(points, 1, 2, 5) == false
    @test Delaunay.is_strictly_right_of(points, 1, 2, 6) == true
    @test Delaunay.is_strictly_right_of(points, 1, 2, 7) == false
end

@testset "Test of Triangle constructors" begin
    points = [
         0  0
         1  0
         1  1
    ]
    @testset "For two vertices" begin
        t = Delaunay.Triangle(1, 2)
        @test t.v3 == 0
        @test t.t1 == t.t2 == t.t3
        next = t.t3
        @test next.v3 == 0
        @test next.v1 == t.v2
        @test next.v2 == t.v1
        @test next.t1 == next.t2 == next.t3 == t
    end
    @testset "For three vertices" begin
        t = Delaunay.Triangle(points, 1, 2, 3)
        # v1 -> v2 -> v3 should be clockwise
        @test t.v1 == 1
        @test t.v2 == 3
        @test t.v3 == 2
        @test Delaunay.is_left_of(points, t.v1, t.v3, t.v2)
        @test t.v1 == 1

        @test t.t1.t3 == t.t2.t3 == t.t3.t3 == t
        @test t.t1.t1 == t.t3.t2 == t.t2
        @test t.t1.t2 == t.t2.t1 == t.t3
        @test t.t2.t2 == t.t3.t1 == t.t1

        @test t.t2.v2 == t.t3.v1 == t.v1
        @test t.t1.v1 == t.t3.v2 == t.v2
        @test t.t1.v2 == t.t2.v1 == t.v3
    end
end
@testset "Test is_hull" begin
    points = [
         0  0
         1  0
         1  1
    ]
    t = Delaunay.Triangle(points, 1, 2, 3)
    @test Delaunay.is_hull(t.t1)
    @test !Delaunay.is_hull(t)
end
@testset "Test hull navigation function" begin
    points = [
         0  0
         1  0
         1  1
    ]
    t = Delaunay.Triangle(points, 1, 2, 3)
    @test Delaunay.next(t.t1) == t.t2
    @test Delaunay.prev(t.t1) == t.t3
end
@testset "Test find_base_tangente" begin
    points = [
         0  0
         1  0
         1  1
         2  0
         3  0
         3  1
    ]
    t1 = Delaunay.Triangle(points, 1, 2, 3)
    t2 = Delaunay.Triangle(points, 4, 5, 6)
    t3 = Delaunay.Triangle(2, 3)
    t4 = Delaunay.Triangle(4, 5)

    t_left, t_right = Delaunay.find_base_tangente(points, t1.t1, t2.t2)
    @test (t_left.v1, t_right.v1) == (2, 4)
end
@testset "Test _is_above_base" begin
    @testset "Test _is_above_base_left" begin
        @testset "For two edges triangles" begin
            points = [
                0 0
                0 1
                1 0
            ]
            t = Delaunay.Triangle(1, 2)
            @test Delaunay._is_above_base_left(t, points, 1, 3) == true
            points = [
                0 -1
                0 0
                1 0
            ]
            t = Delaunay.Triangle(2, 1)
            @test Delaunay._is_above_base_left(t, points, 2, 3) == false
        end
        @testset "For classic triangles" begin
            points = [
                0 0
                1 2
                2 1
                3 0
            ]
            t = Delaunay.Triangle(points, 1, 2, 3)
            @test Delaunay._is_above_base_left(t, points, 1, 4) == true
            points = [
                0 0
                1 -1
                1 -2
                3 0
            ]
            t = Delaunay.Triangle(points, 1, 2, 3)
            @test Delaunay._is_above_base_left(t, points, 1, 4) == false
            points = [
                0 1
                1 0
                2 0
                3 0
            ]
            t = Delaunay.Triangle(points, 1, 2, 3)
            @test Delaunay._is_above_base_left(t, points, 3, 4) == true
            points = [
                0 -1
                1 0
                2 0
                3 0
            ]
            t = Delaunay.Triangle(points, 1, 2, 3)
            @test Delaunay._is_above_base_left(t, points, 3, 4) == false
            points = [
                0 1
                1 1
                2 0
                3 0
            ]
            t = Delaunay.Triangle(points, 1, 2, 3)
            @test Delaunay._is_above_base_left(t, points, 3, 4) == true
            points = [
                0 -1
                1 1
                2 0
                3 0
            ]
            t = Delaunay.Triangle(points, 1, 2, 3)
            @test Delaunay._is_above_base_left(t, points, 3, 4) == false
        end
    end
    @testset "Test _is_above_base_right" begin
        @testset "For two edges triangles" begin
            points = [
                0 0
                1 0
                1 1
            ]
            t = Delaunay.Triangle(2, 3)
            @test Delaunay._is_above_base_right(t, points, 1, 2) == true
            points = [
                0 0
                1 0
                1 -1
            ]
            t = Delaunay.Triangle(3, 2)
            @test Delaunay._is_above_base_right(t, points, 1, 2) == false
        end
        @testset "For classic triangles" begin
            points = [
                0 0
                1 3
                2 1
                3 0
            ]
            t = Delaunay.Triangle(points, 2, 3, 4)
            @test Delaunay._is_above_base_right(t, points, 1, 4) == true
            points = [
                0 0
                1 -3
                2 -2
                3 0
            ]
            t = Delaunay.Triangle(points, 2, 3, 4)
            @test Delaunay._is_above_base_right(t, points, 1, 4) == false
            points = [
                0 0
                1 1
                2 0
                3 1
            ]
            t = Delaunay.Triangle(points, 2, 3, 4)
            @test Delaunay._is_above_base_right(t, points, 1, 3) == true
            points = [
                0 0
                1 -1
                2 0
                3 -1
            ]
            t = Delaunay.Triangle(points, 2, 3, 4)
            @test Delaunay._is_above_base_right(t, points, 1, 3) == false
            points = [
                0 0
                1 0
                2 1
                3 1
            ]
            t = Delaunay.Triangle(points, 2, 3, 4)
            @test Delaunay._is_above_base_right(t, points, 1, 2) == true
            points = [
                0 0
                1 0
                2 -1
                3 -1
            ]
            t = Delaunay.Triangle(points, 2, 3, 4)
            @test Delaunay._is_above_base_right(t, points, 1, 2) == false
        end
    end
end
@testset "Test _next" begin
    @testset "For two edges triangles" begin
        t = Delaunay.Triangle(1, 2)
        @test Delaunay._next_left(t, 1) == t.t3
        @test Delaunay._next_left(Delaunay._next_left(t, 1), 1) == t
        @test Delaunay._next_right(t, 1) == t.t3
        @test Delaunay._next_right(Delaunay._next_right(t, 1), 1) == t
    end
    @testset "For classic triangles" begin
        points = [
            0 0
            1 0
            1 1
        ]
        t = Delaunay.Triangle(points, 1, 2, 3)
        @test Delaunay._next_left(t, 1) == t.t3
        @test Delaunay._next_left(t, 3) == t.t1
        @test Delaunay._next_left(t, 2) == t.t2
        @test Delaunay._next_right(t, 1) == t.t2
        @test Delaunay._next_right(t, 3) == t.t3
        @test Delaunay._next_right(t, 2) == t.t1
    end
end

@testset "I tested every corner case I could think of, but just in case, let's generate 100 graphs and test if they are delaunay 🙈" begin
    for _ in 1:100
        m = Map(20, 20, 10)
        delaunay = Delaunay.is_delaunay(m)
        @test delaunay == true
        if !delaunay
            break
        end
    end
end
end
