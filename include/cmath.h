#pragma once
#include "vectors.h"
#ifdef __cplusplus
extern "C"
{
#endif

float math_std(float *data, int len);
float math_rms(float *data, int len);

/// @brief Interpolate data to a new size.
/// @param originalData The original data.
/// @param originalSize The size of the original data.
/// @param newData A pointer to the new data array.
/// @param newSize The size of the new data array.
/// @return True if successful, false otherwise.
bool math_lin_interpolate(const float *originalData, int originalSize, float *newData, int newSize);
float math_correlate(const float *data1, const float *data2, int size);
void math_normalize(float *data, int size);
void math_zscore(float *data, int size);

/// @brief Scale data to a specific range.
/// @param data The data to scale.
/// @param size Size of the data.
/// @param min The minimum value of the new range.
/// @param max The maximum value of the new range.
void math_scale_to_range(float *data, int size, float min, float max);
/// @brief Remove the DC component from a signal.
/// @param data The data to remove the DC component from.
/// @param size Size of the data.
void math_remove_dc(float *data, int size);
/// @brief Calculate the slope of a line through a given set of points using the least squares method.
/// @param x The x value array
/// @param y The y value array
/// @param n Number of elements in x and y
/// @return The slope of the line
float math_slope(float *x, float *y, int n);
float math_vslope(vector2_t *v, int n);
bool math_zero_cross_x(vector2_t *a, vector2_t *b, float *x);



#ifdef __cplusplus
}
#endif