#pragma once
#include <./glm/glm.hpp>
#include <vector>

struct Particle
{
    Particle()
        : m_Position(0)
        , m_Velocity(0)
        , m_Acceleration(0)
        , m_Color(0)
    {}

    glm::vec3   m_Position; // Center point of particle
    glm::vec3   m_Velocity; // Current particle velocity
    glm::vec3   m_Acceleration; // Current particle acceleration
    glm::vec4   m_Color;    // Particle color
};

typedef std::vector<Particle> ParticleBuffer;