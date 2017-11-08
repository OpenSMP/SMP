//
// Created by riku on 6/15/17.
//

#ifndef CRYPTCONV_CONFIG_H
#define CRYPTCONV_CONFIG_H

#ifdef USE_PARALLEL
#define CRYPTCNN_FOR_PARALLEL _Pragma("omp parallel for")
#else
#define CRYPTCNN_FOR_PARALLEL
#endif
#endif //CRYPTCONV_CONFIG_H
