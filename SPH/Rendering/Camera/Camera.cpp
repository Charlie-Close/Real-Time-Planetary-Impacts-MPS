//
//  Camera.cpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//

#include "Camera.hpp"
#include <iostream>
#include "MatrixMath.h"

Camera::Camera() {
    _perspective = math::makePerspective(45.f * M_PI / 180.f, 1.f, near_plane, far_plane);
    _view = math::makeLookAt(_position, _forward, _up);
}

simd::float4x4 Camera::getMatrix() {
    move();
    _view = math::makeLookAt(_position, _forward, _up);
    return _perspective * _view;
}

void Camera::move() {
    float cameraSpeed = 4.f;
    if (downKeys[0]) {
        _position += cameraSpeed * _forward;
    }
    if (downKeys[1]) {
        _position -= cameraSpeed * _forward;
    }
    if (downKeys[2]) {
        _position -= cameraSpeed * simd::normalize(simd::cross(_forward, _up));
    }
    if (downKeys[3]) {
        _position += cameraSpeed * simd::normalize(simd::cross(_forward, _up));
    }
    if (downKeys[4]) {
        _position += cameraSpeed * _up;
    }
    if (downKeys[5]) {
        _position -= cameraSpeed * _up;
    }
}

void Camera::handleKeyDown(int key) {
    switch (key) {
        case 13: // W
            downKeys[0] = true;
            break;
        case 1: // S
            downKeys[1] = true;
            break;
        case 0: // A
            downKeys[2] = true;
            break;
        case 2: // D
            downKeys[3] = true;
            break;
        case 12: // Q
            downKeys[4] = true;
            break;
        case 14: // E
            downKeys[5] = true;
            break;
    }
}

void Camera::handleKeyUp(int key) {
    switch (key) {
        case 13: // W
            downKeys[0] = false;
            break;
        case 1: // S
            downKeys[1] = false;
            break;
        case 0: // A
            downKeys[2] = false;
            break;
        case 2: // D
            downKeys[3] = false;
            break;
        case 12: // Q
            downKeys[4] = false;
            break;
        case 14: // E
            downKeys[5] = false;
            break;
    }
}

void Camera::handleMouseDown(float x, float y) {
    lastMousePos = { x, y };
}

void Camera::handleMouseDrag(float x, float y) {
    simd_float2 newMousePos = (simd_float2){ x, y };
    simd_float2 offset = newMousePos - lastMousePos;
    lastMousePos = newMousePos;


    // Sensitivity
    float sensitivity = 0.1f;
    yaw += offset.x * sensitivity;
    pitch += offset.y * sensitivity * sensitivity;

    // Constrain pitch
    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    // Update front vector
    simd_float3 direction;
    float yawRad = yaw * M_PI / 180.f;
    direction.x = cos(yawRad) * cos(pitch);
    direction.y = sin(pitch);
    direction.z = sin(yawRad) * cos(pitch);
    _forward = simd::normalize(direction);
}
