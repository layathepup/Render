// takes a 3D mesh and a perspective
// and algebraically transform it into a 2D mesh
// copying the mesh

    // Grammar:
// Mesh         iterable, copyable container of vertices, of edges
    // 2D           vertices, edges span 2 Dimensions
    // 3D           vertices, edges span 3 Dimensions
// Edge         pair of vertex indexes
// Vertex       vector of dimensions
// Dimension    number
// Matrix       vector of vectors
// Perspective  Position vertex, Orientation Vertex, Rotation vector
// Precision    number type, e.g. int, float

#ifndef RENDER_RENDERER_HPP
#define RENDER_RENDERER_HPP

#include <vector>
#include <memory>
#include <cmath>


enum class Dim {
    x, y, z
};

template <typename Precision>
using Vertex = std::vector<Precision>;

// perspective requires a position, orientation, and rotation
template <typename Precision>
using Perspective = std::vector<std::vector<Precision>>;

enum class Perspec {
    Pos, Orient, Rotat
};

template <typename Precision>
using Edge = std::pair<typename std::vector<Precision>::iterator,
                       typename std::vector<Precision>::iterator>;       // indices of mesh

template <typename Precision,
template <typename T, typename A = std::allocator<T>> class Container = std::vector>
struct Mesh {
    Container<Vertex<Precision>> vertices;
    Container<Edge<Precision>> edges;
};


template <typename Elem>
class Matrix_2D {
    // matrix is y-axis vector of x-axis vectors of Elements
    /* ex.  { 1, 0, 0,
              0, 1, 0,
              0, 0, 1 } */
public:
    Matrix_2D(int x, int y)
        :x_size{ x }, y_size{ y }
    {
        //initialize to 0
        for (int i{ 0 }; i<y_size; i++)
            matrix.push_back(std::vector<Elem>(x_size));
    }

    Matrix_2D operator=(std::initializer_list<Elem> lst) {
        if (lst.size() != x_size * y_size)
            throw std::runtime_error("initializer list size does not match matrix size.\n");
        for (Elem elem : lst)
            matrix.push_back(elem);
    }

    std::vector<Elem>& operator[](int i) {
        return matrix[i];
    }

private:
    int x_size;
    int y_size;
    std::vector<std::vector<Elem>> matrix;
};

template <typename Number>
Matrix_2D<Number> create_3x3_rotation_matrix(const Number& theta)
{
    Matrix_2D<Number> mat{3, 3};        // special 3 by 3 matrix
    mat = { std::cos(theta), -std::sin(theta),  0,
            std::sin(theta), std::cos(theta),   0,
            0,               0,                 0   };
    return mat;
}

template <typename Precision>
Precision vector_length(std::vector<Precision> vec)
{
    Precision sum_of_squares{ 0 };
    for (int i{ 0 }; i<vec.size(); i++)
        sum_of_squares += std::pow(vec[i], 2);
    return std::sqrt(sum_of_squares);
}

template <typename Number>
void transform(std::vector<Number>& vec, const Matrix_2D<Number>& mat)
{
    if (vec.size() != mat.x_size() ||
        vec.size() != mat.y_size())
        throw std::runtime_error("cannot transform vector with matrix of different size.\n");

    for (int i{ 0 }; i<vec.size(); i++)
        for (int j{ 0 }; j<vec.size(); j++)
            vec[i] += vec[i] *= mat[i][j];
    return;
}

template <typename Precision,
template <typename T, typename Alloc = std::allocator<T>> class Vertex_Container = std::vector,
template <typename T, typename Alloc = std::allocator<T>> class Mesh_Container = std::vector>
void render(const std::shared_ptr<Mesh_Container<Mesh<Precision, Vertex_Container>>> meshes_ptr,
            const std::shared_ptr<Perspective<Precision>> perspective_ptr)
{
    // Position match
    // translate to perspective is origin
    for (Mesh<Precision, Mesh_Container> mesh : *meshes_ptr)
        for (Vertex<Precision> vertex : mesh.vertices)
            for (int i{ 0 }; i<vertex.size(); i++)
                vertex[i] -= *perspective_ptr[Perspec::Pos][i];

    // Orientation match
    // Orient so that we are looking down the z-axis (depth)
    {
        // calculate angle of xz (ground) projection

        // find length of hypotenuse
        Precision hypo_len_xz = std::sqrt(
        std::exp(*perspective_ptr[Perspec::Orient][Dim::x], 2) +
        std::exp(*perspective_ptr[Perspec::Orient][Dim::z], 2)
        );

        // normalize x and take inverse trig function to find theta
        Precision x = *perspective_ptr[Perspec::Orient][Dim::x] / hypo_len_xz;
        double theta_xz = - (0.5 * M_PI - std::asin(x));

        // create transformation matrix and transform meshes
        for (Mesh<Precision, Mesh_Container> mesh : *meshes_ptr)
            for (Vertex<Precision> vertex : mesh.vertices)
                transform(vertex, create_3x3_rotation_matrix(theta_xz));

        // same process for yz projection (does this repetition justify its own function?)
        Precision hypo_len_yz = std::sqrt(
                std::exp(*perspective_ptr[Perspec::Orient][Dim::y], 2) +
                std::exp(*perspective_ptr[Perspec::Orient][Dim::z], 2)
        );

        Precision y = *perspective_ptr[Perspec::Orient][Dim::y] / hypo_len_yz;
        double theta_yz = -std::asin(y);

        for (Mesh<Precision, Mesh_Container> mesh : *meshes_ptr)
            for (Vertex<Precision> vertex : mesh.vertices)
                transform(vertex, create_3x3_rotation_matrix(theta_xz));
    }

    // Rotation match
    // Rotate so that the rotation perspective vector is (0, 0, 1)
    {
        Precision hypo_len_xy = std::sqrt(
                std::exp(*perspective_ptr[Perspec::Rotat][Dim::x], 2) +
                std::exp(*perspective_ptr[Perspec::Rotat][Dim::y], 2)
        );

        Precision x = *perspective_ptr[Perspec::Rotat][Dim::x] / hypo_len_xy;
        double theta_xy = - (0.5 * M_PI - std::asin(x));

        for (Mesh<Precision, Mesh_Container> mesh : *meshes_ptr)
            for (Vertex<Precision> vertex : mesh.vertices)
                transform(vertex, create_3x3_rotation_matrix(theta_xy));
    }

}



template <typename Precision,
template <typename T, typename Alloc = std::allocator<T>> class Container = std::vector>
class Renderer {
public:
    Renderer(const std::shared_ptr<Mesh<Precision, Container>> const mesh_ptr,
             const Vertex<Precision>& perspective)
        :mesh_ptr{ mesh_ptr } {};



private:
    const std::shared_ptr<Mesh<Precision, Container>> const mesh_ptr;
};


#endif //RENDER_RENDERER_HPP
