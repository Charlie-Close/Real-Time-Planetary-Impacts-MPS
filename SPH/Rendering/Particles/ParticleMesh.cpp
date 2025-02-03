//
//  ParticleMesh.cpp
//  SPH
//
//  Created by Charlie Close on 24/01/2025.
//

#include "ParticleMesh.hpp"
#include <unordered_map>

struct Vec3Hash {
    std::size_t operator()(const simd_float3& v) const {
        return std::hash<float>()(v.x) ^ std::hash<float>()(v.y) ^ std::hash<float>()(v.z);
    }
};

// Custom equality operator for simd_float3
struct Vec3Equal {
    bool operator()(const simd_float3& a, const simd_float3& b) const {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }
};

void reduce_mesh_data(std::vector<simd_float3>& verts,
                      std::vector<std::vector<uint16_t>>& faces,
                      std::vector<simd_float3>& normals) {
    
     std::unordered_map<simd_float3, int, Vec3Hash, Vec3Equal>  unique_verts_map;
     std::vector<simd_float3> unique_verts;
     std::vector<int> inverse_indices(verts.size());

    // Find unique vertices and map to new indices
    for (size_t i = 0; i < verts.size(); ++i) {
        const simd_float3& vert = verts[i];
        if (unique_verts_map.find(vert) == unique_verts_map.end()) {
            unique_verts_map[vert] = unique_verts.size();
            unique_verts.push_back(vert);
        }
        inverse_indices[i] = unique_verts_map[vert];
    }

    // Map faces to new indices from unique vertices
    std::vector<std::vector<uint16_t>> new_faces(faces.size(), std::vector<uint16_t>(3));
    for (size_t i = 0; i < faces.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            new_faces[i][j] = inverse_indices[faces[i][j]];
        }
    }

    // Accumulate data for unique vertices
    std::vector<simd_float3> normal_accumulator(unique_verts.size(), simd_float3(0.0f));

    for (size_t i = 0; i < normals.size(); ++i) {
        normal_accumulator[inverse_indices[i]] += normals[i];
    }

    // Normalize the accumulated normals
    std::vector<simd_float3> avg_normals(unique_verts.size());

    for (size_t i = 0; i < normal_accumulator.size(); ++i) {
        float length = simd::length(normal_accumulator[i]);
        avg_normals[i] = (length > 1e-8) ? (normal_accumulator[i] / length) : simd_float3(0.0f);
    }

    verts = unique_verts;
    faces = new_faces;
    normals = avg_normals;
}

// Function to compute normals for vertices in a mesh
std::vector<simd_float3> compute_normals(const std::vector<simd_float3>& vertices,
                                       const std::vector<std::vector<uint16_t>>& faces,
                                       const simd_float3 center) {
    std::vector<simd_float3> normals(vertices.size(), simd_float3(0.0f));

    // Iterate over all faces
    for (const auto& face : faces) {
        int i0 = face[0];
        int i1 = face[1];
        int i2 = face[2];

        // Compute the two edge vectors of the triangle
        simd_float3 v0 = vertices[i1] - vertices[i0];
        simd_float3 v1 = vertices[i2] - vertices[i0];

        // Compute the normal of the triangle (cross product of edges)
        simd_float3 normal = simd::cross(v0, v1);

        // If a center is provided, ensure normal points toward it
        simd_float3 view = center - vertices[i0];
        if (simd::dot(normal, view) < 0) {
            normal = -normal;
        }

        // Normalize the normal vector
        float length = simd::length(normal);
        if (length > 0) {
            normal /= length;
        }

        // Add this normal to each vertex of the triangle
        normals[i0] += normal;
        normals[i1] += normal;
        normals[i2] += normal;
    }

    // Normalize all vertex normals
    for (auto& normal : normals) {
        float length = simd::length(normal);
        if (length > 0) {
            normal /= length;
        }
    }

    return normals;
}

// Helper function to subdivide each triangle face in the mesh
void subdivide(std::vector<simd_float3>& vertices, std::vector<std::vector<uint16_t>>& faces) {
    std::vector<std::vector<uint16_t>> new_faces;
    std::unordered_map<uint64_t, uint16_t> midpoint_cache;

    auto midpoint = [&](uint16_t v1, uint16_t v2) -> int {
        uint64_t key = (static_cast<uint64_t>(std::min(v1, v2)) << 32) | std::max(v1, v2);
        if (midpoint_cache.count(key)) return midpoint_cache[key];

        simd_float3 mid = simd::normalize((vertices[v1] + vertices[v2]) * 0.5f);
        vertices.push_back(mid);
        return midpoint_cache[key] = vertices.size() - 1;
    };

    for (const auto& face : faces) {
        uint16_t v1 = face[0], v2 = face[1], v3 = face[2];
        uint16_t a = midpoint(v1, v2);
        uint16_t b = midpoint(v2, v3);
        uint16_t c = midpoint(v3, v1);

        new_faces.push_back({v1, a, c});
        new_faces.push_back({v2, b, a});
        new_faces.push_back({v3, c, b});
        new_faces.push_back({a, b, c});
    }

    faces = new_faces;
}


// Helper function to create an initial icosphere, which can then be subdivided
std::tuple<std::vector<simd_float3>, std::vector<std::vector<uint16_t>>> generate_icosphere(int subdivisions) {
    const float X = 0.525731f;
    const float Z = 0.850651f;

    std::vector<simd_float3> vertices = {
        (simd_float3){-X, 0.0f, Z}, (simd_float3){X, 0.0f, Z}, (simd_float3){-X, 0.0f, -Z}, (simd_float3){X, 0.0f, -Z},
        (simd_float3){0.0f, Z, X}, (simd_float3){0.0f, Z, -X}, (simd_float3){0.0f, -Z, X}, (simd_float3){0.0f, -Z, -X},
        (simd_float3){Z, X, 0.0f}, (simd_float3){-Z, X, 0.0f}, (simd_float3){Z, -X, 0.0f}, (simd_float3){-Z, -X, 0.0f}
    };

    std::vector<std::vector<uint16_t>> faces = {
        {0, 4, 1}, {0, 9, 4}, {9, 5, 4}, {4, 5, 8}, {4, 8, 1},
        {8, 10, 1}, {8, 3, 10}, {5, 3, 8}, {5, 2, 3}, {2, 7, 3},
        {7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6},
        {6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5}, {7, 2, 11}
    };

    // Subdivide faces to increase the number of vertices
    for (int i = 0; i < subdivisions; ++i) {
        subdivide(vertices, faces);
    }

    return {vertices, faces};
}


std::tuple<std::vector<simd_float3>, std::vector<simd_float3>, std::vector<uint16_t>> generateSphere(float size, int subdivisions) {
    std::vector<simd_float3> vertices;
    std::vector<std::vector<uint16_t>> faces;
    std::tie(vertices, faces) = generate_icosphere(subdivisions);
    std::vector<simd_float3> normals = compute_normals(vertices, faces, (simd_float3){ 0, 0, 0 });
    reduce_mesh_data(vertices, faces, normals);
    
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] *= size;
    }
    
    std::vector<uint16_t> indices(faces.size() * 3);
    for (int i = 0; i < faces.size(); i ++) {
        indices[i * 3] = faces[i][0];
        indices[i * 3 + 1] = faces[i][1];
        indices[i * 3 + 2] = faces[i][2];
    }

    return {vertices, normals, indices};
}
