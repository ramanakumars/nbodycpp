#pragma once

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <omp.h>

#define CHUNK_SIZE 500
const double M_PI = 3.14159;

const double GRAV_G = 1;
const double MASS_REF = 0.1;
const double ALPHA = 0.5;

enum transport_type {
    LEAPFROG,
    YOSHIDA
};

extern transport_type TRANSPORT_TYPE;
