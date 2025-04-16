#include <stdbool.h>
#include <string.h>
#include <float.h> // DBL_EPSILON
#include "math.h"
#include "cmath.h"
#include "vectors.h"

extern "C"
{

    /// @brief Calculate mean of data.
    /// mean = 1/N * sum(x)
    /// @param data 
    /// @param len 
    /// @return 
    float math_mean(const float* data, int len)
    {
        float mean = 0.0f;
        for (int i = 0; i < len; i++)
            mean += data[i];
        mean /= len;
        return mean;
    }

    /// @brief Calculate standard deviation of data.
    /// std = sqrt(1/N * sum((x - mean)^2))
    /// @param data 
    /// @param len 
    /// @return 
    float math_std(float* data, int len)
    {
        float std = 0.0f;
        float mean = math_mean(data, len);

        // Calculate distance from mean, square and sum
        for (int i = 0; i < len; i++)
        {
            std += (data[i] - mean) * (data[i] - mean);
        }
        std /= len;

        // Square root
        std = sqrtf(std);
        return std;
    }

    float math_rms(float* data, int len)
    {
        float rms = 0.0f;
        for (int i = 0; i < len; i++)
        {
            rms += data[i] * data[i];
        }
        rms /= len;
        rms = sqrtf(rms);
        return rms;
    }

    bool math_lin_interpolate(const float* originalData, int originalSize, float* newData, int newSize)
    {
        if (originalSize < 2 || newSize < 2) {
            return false;
        }
        if (originalSize == newSize) {
            memcpy(newData, originalData, originalSize * sizeof(float));
            return true;
        }
        for (int i = 0; i < newSize; ++i) {
            double pos = ((double)i * (originalSize - 1)) / (newSize - 1);
            int posInt = (int)pos;
            double frac = pos - posInt;

            if (posInt >= originalSize - 1) {
                newData[i] = originalData[originalSize - 1];
            }
            else {
                newData[i] = (1 - frac) * originalData[posInt] + frac * originalData[posInt + 1];
            }
        }
        return true;
    }

    float math_correlate(const float* data1, const float* data2, int size)
    {
        if (size < 2) {
            // Not enough data to compute correlation
            return 0.0f;
        }

        double mean1 = 0, mean2 = 0;
        double sum1 = 0, sum2 = 0, sumProduct = 0;
        int i;

        // Calculate means of the two datasets
        for (i = 0; i < size; i++) {
            mean1 += data1[i];
            mean2 += data2[i];
        }
        mean1 /= size;
        mean2 /= size;

        // Calculate standard deviations and covariance
        for (i = 0; i < size; i++) {
            double delta1 = data1[i] - mean1;
            double delta2 = data2[i] - mean2;
            sum1 += delta1 * delta1;
            sum2 += delta2 * delta2;
            sumProduct += delta1 * delta2;
        }

        // Use sample standard deviation (n - 1)
        double stdev1 = sqrt(sum1 / (size - 1));
        double stdev2 = sqrt(sum2 / (size - 1));

        // Protect against division by zero
        if (stdev1 < DBL_EPSILON || stdev2 < DBL_EPSILON) {
            return 0.0f;
        }

        return (float)(sumProduct / ((size - 1) * stdev1 * stdev2));
    }

    void math_normalize(float* data, int size)
    {
        float min = data[0];
        float max = data[0];
        for (int i = 1; i < size; i++) {
            if (data[i] < min) {
                min = data[i];
            }
            if (data[i] > max) {
                max = data[i];
            }
        }
        float range = max - min;
        for (int i = 0; i < size; i++) {
            data[i] = (data[i] - min) / range;
        }
    }


    void math_zscore(float* data, int size)
    {
        float mean = math_mean(data, size);
        float std = math_std(data, size);
        for (int i = 0; i < size; i++) {
            data[i] = (data[i] - mean) / std;
        }
    }


    void math_remove_dc(float* data, int size)
    {
        float mean = math_mean(data, size);
        for (int i = 0; i < size; i++) {
            data[i] -= mean;
        }
    }


    void math_scale_to_range(float* data, int size, float min, float max)
    {
        float dataMin = data[0];
        float dataMax = data[0];
        for (int i = 1; i < size; i++) {
            if (data[i] < dataMin) {
                dataMin = data[i];
            }
            if (data[i] > dataMax) {
                dataMax = data[i];
            }
        }
        float dataRange = dataMax - dataMin;
        float scale = (max - min) / dataRange;
        for (int i = 0; i < size; i++) {
            data[i] = (data[i] - dataMin) * scale + min;
        }
    }

    float math_slope(float* x, float* y, int n)
    {
        float sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x_squared = 0.0;

        // Calculate the sums needed for the slope formula
        for (int i = 0; i < n; i++) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            sum_x_squared += x[i] * x[i];
        }

        // Apply the formula for slope
        float slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x * sum_x);

        return slope;
    }

    float math_vslope(vector2_t* v, int n)
    {
        float sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x_squared = 0.0;

        // Calculate the sums needed for the slope formula
        for (int i = 0; i < n; i++) {
            vector2_t* v2 = &v[i];
            sum_x += v2->x;
            sum_y += v2->y;
            sum_xy += v2->x * v2->y;
            sum_x_squared += v2->x * v2->x;
        }

        // Apply the formula for slope
        float slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x * sum_x);

        return slope;
    }
} // extern "C"


