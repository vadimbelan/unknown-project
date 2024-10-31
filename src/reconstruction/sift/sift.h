namespace reconstruction {
    class SIFT {
    public:
        SIFT(double contrast_threshold = 0.5) : contrast_threshold(contrast_threshold) {}
        void detectAndCompute();
    protected:
        void buildPyramids();
        void findLocalExtremasAndDescribe();
        bool buildLocalOrientationHists();
        bool buildDescriptor();
        double contrast_threshold;
    };
}
