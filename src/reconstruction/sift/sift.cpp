#include "sift.h"
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <libutils/rasserts.h>

#define NUM_OCTAVES 3
#define LAYERS_PER_OCTAVE 3
#define GAUSSIAN_IMAGES_PER_OCTAVE (LAYERS_PER_OCTAVE + 3)
#define DOG_IMAGES_PER_OCTAVE (LAYERS_PER_OCTAVE + 2)
#define INITIAL_SIGMA 0.75
#define INPUT_BLUR_SIGMA 1.0
#define NUM_ORIENTATION_HIST_BINS 36
#define ORIENTATION_WINDOW_RADIUS 3
#define ORIENTATION_PEAK_RATIO 0.80
#define DESCRIPTOR_GRID_SIZE 4
#define DESCRIPTOR_HIST_BINS 8
#define DESCRIPTOR_SAMPLE_COUNT 4
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

void reconstruction::SIFT::findLocalExtremasAndDescribe(const std::vector<cv::Mat> &gaussianPyramid, const std::vector<cv::Mat> &DoGPyramid,
                                             std::vector<cv::KeyPoint> &keyPoints, cv::Mat &desc) {
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

bool reconstruction::SIFT::buildLocalOrientationHists(const cv::Mat &img, size_t i, size_t j, size_t radius,
                                           std::vector<float> &votes, float &biggestVote) {
    votes.resize(NUM_ORIENTATION_HIST_BINS, 0.0f);
    biggestVote = 0.0;

    if (i < radius || i + radius >= img.cols || j < radius || j + radius >= img.rows)
        return false;

    float sum[NUM_ORIENTATION_HIST_BINS] = {0.0f};

    for (size_t y = j - radius; y <= j + radius; ++y) {
        for (size_t x = i - radius; x <= i + radius; ++x) {
            float dx = img.at<float>(y, x + 1) - img.at<float>(y, x - 1);
            float dy = img.at<float>(y + 1, x) - img.at<float>(y - 1, x);

            double magnitude = sqrt(dx * dx + dy * dy);
            double orientation = atan2(dy, dx) * 180.0 / M_PI;
            orientation = fmod(orientation + 360.0, 360.0);

            size_t bin = static_cast<size_t>(orientation / (360.0 / NUM_ORIENTATION_HIST_BINS));
            rassert(bin < NUM_ORIENTATION_HIST_BINS, 361236315613);

            sum[bin] += magnitude;
        }
    }

    for (size_t bin = 0; bin < NUM_ORIENTATION_HIST_BINS; ++bin) {
        float prev = sum[(bin + NUM_ORIENTATION_HIST_BINS - 1) % NUM_ORIENTATION_HIST_BINS];
        float next = sum[(bin + 1) % NUM_ORIENTATION_HIST_BINS];
        votes[bin] = (prev + sum[bin] + next) / 3.0f;
        biggestVote = std::max(biggestVote, votes[bin]);
    }

    return true;
}

bool reconstruction::SIFT::buildDescriptor(const cv::Mat &img, float px, float py, double descrSampleRadius, float angle,
                                std::vector<float> &descriptor) {
    cv::Mat relativeShiftRotation = cv::getRotationMatrix2D(cv::Point2f(0.0f, 0.0f), -angle, 1.0);
    const double smpW = 2.0 * descrSampleRadius - 1.0;

    descriptor.resize(DESCRIPTOR_GRID_SIZE * DESCRIPTOR_GRID_SIZE * DESCRIPTOR_HIST_BINS, 0.0f);
    for (int hstj = 0; hstj < DESCRIPTOR_GRID_SIZE; ++hstj) {
        for (int hsti = 0; hsti < DESCRIPTOR_GRID_SIZE; ++hsti) {
            float sum[DESCRIPTOR_HIST_BINS] = {0.0f};

            for (int smpj = 0; smpj < DESCRIPTOR_SAMPLE_COUNT; ++smpj) {
                for (int smpi = 0; smpi < DESCRIPTOR_SAMPLE_COUNT; ++smpi) {
                    cv::Point2f shift(((-DESCRIPTOR_GRID_SIZE / 2.0 + hsti) * DESCRIPTOR_SAMPLE_COUNT + smpi) * smpW,
                                      ((-DESCRIPTOR_GRID_SIZE / 2.0 + hstj) * DESCRIPTOR_SAMPLE_COUNT + smpj) * smpW);
                    std::vector<cv::Point2f> shiftInVector(1, shift);
                    cv::transform(shiftInVector, shiftInVector, relativeShiftRotation);
                    shift = shiftInVector[0];

                    int x = static_cast<int>(px + shift.x);
                    int y = static_cast<int>(py + shift.y);

                    if (y - 1 < 0 || y + 1 >= img.rows || x - 1 < 0 || x + 1 >= img.cols)
                        return false;

                    float dx = img.at<float>(y, x + 1) - img.at<float>(y, x - 1);
                    float dy = img.at<float>(y + 1, x) - img.at<float>(y - 1, x);
                    double magnitude = sqrt(dx * dx + dy * dy);

                    double orientation = atan2(dy, dx) * 180.0 / M_PI - angle;
                    orientation = fmod(orientation + 360.0, 360.0);

                    float binFloat = orientation / (360.0 / DESCRIPTOR_HIST_BINS);
                    int bin = static_cast<int>(binFloat);
                    float weightBin = binFloat - bin;
                    sum[bin] += (1.0f - weightBin) * magnitude;
                    sum[(bin + 1) % DESCRIPTOR_HIST_BINS] += weightBin * magnitude;
                }
            }

            float *votes = &(descriptor[(hstj * DESCRIPTOR_GRID_SIZE + hsti) * DESCRIPTOR_HIST_BINS]);
            for (int bin = 0; bin < DESCRIPTOR_HIST_BINS; ++bin) {
                votes[bin] = sum[bin];
            }
        }
    }

    float norm = 0.0f;
    for (float val : descriptor) norm += val * val;
    norm = sqrt(norm);
    for (float &val : descriptor) val /= norm;

    for (float &val : descriptor) {
        if (val > 0.2f) val = 0.2f;
    }
    norm = 0.0f;
    for (float val : descriptor) norm += val * val;
    norm = sqrt(norm);
    for (float &val : descriptor) val /= norm;

    return true;
}
