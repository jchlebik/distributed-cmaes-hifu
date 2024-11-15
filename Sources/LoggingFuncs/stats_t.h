#pragma once

/**
 * @brief a structure to contain relevant statistics information to create a boxplot. 
 * 
 */
typedef struct stats_t 
{
    double q1;
    double q3;
    double median;
    double mean;
    double irq;
    double min;
    double max;
} stats;
