#include <vector>
#include <opencv2/core.hpp>

namespace reconstruction {
    class SIFT {
    public:
        SIFT(double contrast_threshold = 0.5) : contrast_threshold(contrast_threshold) {}
        void detectAndCompute();
    protected:
        void buildPyramids(const cv::Mat &imgOrg, std::vector<cv::Mat> &gaussianPyramid, std::vector<cv::Mat> &DoGPyramid);
        void findLocalExtremasAndDescribe();
        bool buildLocalOrientationHists();
        bool buildDescriptor();
        double contrast_threshold;
    };

}
