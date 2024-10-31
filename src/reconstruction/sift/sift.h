#include <vector>
#include <opencv2/core.hpp>

namespace reconstruction {
    class SIFT {
    public:
        SIFT(double contrast_threshold = 0.5) : contrast_threshold(contrast_threshold) {}
        void detectAndCompute();
    protected:
        void buildPyramids(const cv::Mat &imgOrg, std::vector<cv::Mat> &gaussianPyramid, std::vector<cv::Mat> &DoGPyramid);
        void findLocalExtremasAndDescribe(const std::vector<cv::Mat> &gaussianPyramid, const std::vector<cv::Mat> &DoGPyramid,
                                          std::vector<cv::KeyPoint> &keyPoints, cv::Mat &desc);
        bool buildLocalOrientationHists(const cv::Mat &img, size_t i, size_t j, size_t radius,
                                        std::vector<float> &votes, float &biggestVote);
        bool buildDescriptor(const cv::Mat &img, float px, float py, double descrRadius, float angle,
                             std::vector<float> &descriptor);
        double contrast_threshold;
    };
}
