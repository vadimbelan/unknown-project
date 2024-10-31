#include "sift.h"
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <libutils/rasserts.h>

#define NUM_OCTAVES 3 +
#define LAYERS_PER_OCTAVE 3 +
#define GAUSSIAN_IMAGES_PER_OCTAVE (LAYERS_PER_OCTAVE + 3) +
#define DOG_IMAGES_PER_OCTAVE (LAYERS_PER_OCTAVE + 2) +
#define INITIAL_SIGMA 0.75 +
#define INPUT_BLUR_SIGMA ?
#define NUM_ORIENTATION_HIST_BINS 36
#define ORIENTATION_WINDOW_RADIUS 3
#define ORIENTATION_PEAK_RATIO 0.80
#define DESCRIPTOR_GRID_SIZE 4
#define DESCRIPTOR_HIST_BINS 8
#define DESCRIPTOR_SAMPLE_COUNT ?
#define DESCRIPTOR_WINDOW_RADIUS 1.0

void reconstruction::ScaleInvariantFeatureTransform::detectAndCompute() { }

void reconstruction::SIFT::buildPyramids(const cv::Mat &imgOrg, std::vector<cv::Mat> &gaussianPyramid, std::vector<cv::Mat> &DoGPyramid) {
    gaussianPyramid.resize(NUM_OCTAVES * GAUSSIAN_IMAGES_PER_OCTAVE);
    const double k = pow(2.0, 1.0 / LAYERS_PER_OCTAVE);

    for (size_t octave = 0; octave < NUM_OCTAVES; ++octave) {
        int layer = 0;
        if (octave == 0) {
            gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE + layer] = imgOrg.clone();
        } else {
            size_t prevOctave = octave - 1;
            cv::Mat img = gaussianPyramid[prevOctave * GAUSSIAN_IMAGES_PER_OCTAVE + (GAUSSIAN_IMAGES_PER_OCTAVE - 3)].clone();
            cv::resize(img, img, cv::Size(), 0.5, 0.5, cv::INTER_NEAREST);
            gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE + layer] = img;
        }

        #pragma omp parallel for
        for (ptrdiff_t layer = 1; layer < GAUSSIAN_IMAGES_PER_OCTAVE; ++layer) {
            double sigmaCur  = INITIAL_SIGMA * pow(2.0, octave) * pow(k, layer);
            double sigmaFirstLayer = INITIAL_SIGMA * pow(2.0, octave);
            double sigma = sqrt(sigmaCur * sigmaCur - sigmaFirstLayer * sigmaFirstLayer);

            cv::Mat imgLayer = gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE].clone();
            cv::GaussianBlur(imgLayer, imgLayer, cv::Size(0, 0), sigma, sigma);
            gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE + layer] = imgLayer;
        }
    }

    DoGPyramid.resize(NUM_OCTAVES * DOG_IMAGES_PER_OCTAVE);

    #pragma omp parallel for
    for (ptrdiff_t octave = 0; octave < NUM_OCTAVES; ++octave) {
        for (size_t layer = 1; layer < GAUSSIAN_IMAGES_PER_OCTAVE; ++layer) {
            int prevLayer = layer - 1;
            cv::Mat imgPrevGaussian = gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE + prevLayer];
            cv::Mat imgCurGaussian  = gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE + layer];
            cv::Mat imgCurDoG = imgCurGaussian - imgPrevGaussian;

            int dogLayer = layer - 1;
            DoGPyramid[octave * DOG_IMAGES_PER_OCTAVE + dogLayer] = imgCurDoG;
        }
    }
}

void reconstruction::SIFT::findLocalExtremasAndDescribe() {
    std::vector<std::vector<float>> pointsDesc;

    #pragma omp parallel
    {
        std::vector<cv::KeyPoint> thread_points;
        std::vector<std::vector<float>> thread_descriptors;

        for (size_t octave = 0; octave < NUM_OCTAVES; ++octave) {
            double octave_downscale = pow(2.0, octave);
            for (size_t layer = 1; layer + 1 < DOG_IMAGES_PER_OCTAVE; ++layer) {
                const cv::Mat prev = DoGPyramid[octave * DOG_IMAGES_PER_OCTAVE + layer - 1];
                const cv::Mat cur  = DoGPyramid[octave * DOG_IMAGES_PER_OCTAVE + layer];
                const cv::Mat next = DoGPyramid[octave * DOG_IMAGES_PER_OCTAVE + layer + 1];
                const cv::Mat DoGs[3] = {prev, cur, next};

                #pragma omp for
                for (ptrdiff_t j = 1; j < cur.rows - 1; ++j) {
                    for (ptrdiff_t i = 1; i < cur.cols - 1; ++i) {
                        float center = DoGs[1].at<float>(j, i);
                        bool is_extremum = true;

                        for (int dz = -1; dz <= 1 && is_extremum; ++dz) {
                            for (int dy = -1; dy <= 1 && is_extremum; ++dy) {
                                for (int dx = -1; dx <= 1 && is_extremum; ++dx) {
                                    if (dz != 0 || dy != 0 || dx != 0) {
                                        float neighbor = DoGs[1 + dz].at<float>(j + dy, i + dx);
                                        if ((neighbor >= center && center > 0) || (neighbor <= center && center < 0)) {
                                            is_extremum = false;
                                        }
                                    }
                                }
                            }
                        }
                        if (!is_extremum) continue;

                        cv::KeyPoint kp;
                        float contrast = center;
                        if (contrast < contrast_threshold) continue;

                        kp.pt = cv::Point2f((i + 0.5) * octave_downscale, (j + 0.5) * octave_downscale);
                        kp.response = fabs(contrast);

                        const double k = pow(2.0, 1.0 / LAYERS_PER_OCTAVE);
                        double sigmaCur = INITIAL_SIGMA * pow(2.0, octave) * pow(k, layer);
                        kp.size = 2.0 * sigmaCur * 5.0;

                        cv::Mat img = gaussianPyramid[octave * GAUSSIAN_IMAGES_PER_OCTAVE + layer];
                        std::vector<float> votes;
                        float biggestVote;
                        int oriRadius = (int) (ORIENTATION_WINDOW_RADIUS * (1.0 + k * (layer - 1)));
                        if (!buildLocalOrientationHists(img, i, j, oriRadius, votes, biggestVote)) continue;

                        for (size_t bin = 0; bin < NUM_ORIENTATION_HIST_BINS; ++bin) {
                            float prevValue = votes[(bin + NUM_ORIENTATION_HIST_BINS - 1) % NUM_ORIENTATION_HIST_BINS];
                            float value = votes[bin];
                            float nextValue = votes[(bin + 1) % NUM_ORIENTATION_HIST_BINS];
                            if (value > prevValue && value > nextValue && votes[bin] > biggestVote * ORIENTATION_PEAK_RATIO) {
                                kp.angle = fmod(bin * (360.0 / NUM_ORIENTATION_HIST_BINS) + 360.0, 360.0);

                                std::vector<float> descriptor;
                                double descrSampleRadius = (DESCRIPTOR_WINDOW_RADIUS * (1.0 + k * (layer - 1)));
                                if (!buildDescriptor(img, kp.pt.x, kp.pt.y, descrSampleRadius, kp.angle, descriptor)) continue;

                                thread_points.push_back(kp);
                                thread_descriptors.push_back(descriptor);
                            }
                        }
                    }
                }
            }
        }

        #pragma omp critical
        {
            keyPoints.insert(keyPoints.end(), thread_points.begin(), thread_points.end());
            pointsDesc.insert(pointsDesc.end(), thread_descriptors.begin(), thread_descriptors.end());
        }
    }

    rassert(pointsDesc.size() == keyPoints.size(), 12356351235124);
    desc = cv::Mat(pointsDesc.size(), DESCRIPTOR_GRID_SIZE * DESCRIPTOR_GRID_SIZE * DESCRIPTOR_HIST_BINS, CV_32FC1);
    for (size_t j = 0; j < pointsDesc.size(); ++j) {
        rassert(pointsDesc[j].size() == DESCRIPTOR_GRID_SIZE * DESCRIPTOR_GRID_SIZE * DESCRIPTOR_HIST_BINS, 1253351412421);
        for (size_t i = 0; i < pointsDesc[j].size(); ++i) {
            desc.at<float>(j, i) = pointsDesc[j][i];
        }
    }
}

bool reconstruction::SIFT::buildDescriptor() { 
// Создание матрицы поворота для учета угла ориентации ключевой точки
// Вычисление ширины окна выборки
// Инициализация вектора для хранения дескриптора

// Цикл по всем ячейкам дескриптора
        // Инициализация массива для суммирования значений гистограммы
        // Цикл по всем выборкам в ячейке дескриптора
                // Вычисление смещения относительно центра ячейки
                // Применение поворота к смещению
                // Вычисление координат пикселя на изображении
                // Проверка, что пиксель находится внутри изображения
                // Вычисление градиента по оси X
                // Вычисление градиента по оси Y           
                // Вычисление величины градиента
                // Вычисление ориентации градиента в градусах с учетом угла ориентации ключевой точки
                // Нормализация ориентации в диапазон [0, 360)
                // Вычисление бина гистограммы для текущей ориентации
                // Вычисление веса для интерполяции между бинами
                // Добавление величины градиента в соответствующий бин гистограммы
        // Запись значений гистограммы в дескриптор

// Нормализация дескриптора
// Ограничение максимального значения в дескрипторе
// Повторная нормализация дескриптора
// Возвращение успешного результата
}
