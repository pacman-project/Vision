#include <unordered_map>
#include <fstream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <iostream>
#include <limits>
#include <bitset>   //for hashing
#include <cmath>
#include <Eigen/Dense>

#include <boost/filesystem.hpp>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>   //for colormap


using namespace cv;
using namespace std;
using namespace boost::filesystem;
class FirstLayer{
public:
    ///Constructor
    FirstLayer();
    ///Destructor
    ~FirstLayer();
    /**
    * @brief Writes a matrix to the text file.
    *
    * This method is storing a matrix in a text file in order to get better
    * view at the exact values at certain locations.
    * @param m matrix to be stored in the text file
    * @param file name of the text file
    * @return 0 for success return
    */
    int matrixToTxtFile(const Mat &m, string filename);
    /**
    * @brief Writes a vector of vectors to the text file.
    *
    * This method is storing a vector of vectors in a text file in order to get better
    * view at the exact values at certain locations.
    * @param data vector to be stored in the text file
    * @param file name of the text file
    * @return 0 for success return
    */
    int vectorToTxtFile(const vector<vector<double> > data, string filename);
    //templating required
    /**
    * @brief Writes a vector of vectors to the text file.
    *
    * This method is storing a vector of vectors in a text file in order to get better
    * view at the exact values at certain locations.
    * @param data vector to be stored in the text file
    * @param file name of the text file
    * @return 0 for success return
    */
    int vectorToTxtFile(const vector<vector<float> > data, string filename);
    /**
    * @brief Calculates the gradient of the depth image.
    *
    * This method is calculating the gradient of the image using SCHARR mask.
    * Gradient orientation is calculated but not its magnitude.
    * Method includes initial smoothing with Gaussian Filter.
    * @param m input matrix (image)
    * @param gradX storage for gradient X-gradient
    * @param gradY storage for gradient Y-gradient
    * @param gradOrient storage for gradient orientation
    * @param gradMagnit storage for gradient magnitude
    * @return 0 for success return
    */
    int calcGradient(const Mat &m, Mat &gradX, Mat &gradY, Mat &gradOrient, Mat &gradMagnit);
    /**
    * @brief Calculates histogram of the depth image.
    *
    * This method is calculating the Histogram of the gradient.
    * All the parameters are set using class variables.
    * @param m input matrix (image)
    * @param hist storage for a histogram
    * @return 0 for success return
    */
    int calcHistogram(const Mat &m, Mat &hist);
    /**
    * @brief Calculates color mapping for displaying depth images.
    *
    * This method is map 16-bit image into color map for visualization purposes
    * @param m input matrix (image)
    * @param colorMap storage for a histogram
    * @return 0 for success return
    */
    int calcColorMap(const Mat &m, Mat &colorMap);

    /**
    * @brief Calculates receptive field without overlap.
    *
    * Receptive field at the botom layer is calculated as a mean of values in
    * the neighbourhood specified by the receptive field size.
    * Border is replicated on the right and bottom part of the image.
    * @param m input matrix (image)
    * @param recept storage for an image after contraction
    * @param recFieldSize size of the receptive field
    * @return 0 for success return
    */
    int calcRecptField(const Mat &m, Mat &recept, int recFieldSize);
    /**
    * @brief Segments parts according to distance measure
    *
    * Use this function together with calcDistanceMeasure.
    * Segmentation is performed according to distance measure defined in calcDistanceMeasure.
    * List of parts obtained from calcDistanceMeasure is stored in
    * _parts[layerNumber][partNumber] temporary to count all the parts
    * Store the relation on the current layer _layerMap[layerNumber][partNumber];
    * @param data calculated distance measure for all the parts
    * @param output storage for the image after segmentation
    * @param layerNum number of layer which is currently processed
    * @param difference threshold used for comparing similarity of parts
    * @return 0 for success return
    */
    int segmentParts(vector<vector<double> > &data, Mat &output, int layerNum, int difference);


    /**
    * @brief Calculates distance between patches at the first layer
    *
    * This method allows to obtain the distances between all the patches of size
    * (receptiveField x receptiveField)present in the image of the first layer.
    * Only at the first layer calculation are performed on the on actual data.
    * @param m input matrix (image of the selected layer)
    * @param data storage for calculated distance measure for all the parts (output)
    * @param recFieldSize size of the receptive field for the current layer
    * @return 0 for success return
    */
    int calcDistanceMeasureHalf(const Mat &m, vector<vector<double> > &data, int recFieldSize);

    /**
    * @brief Calculates distance between patches at the first layer
    *
    * This method allows to obtain the distances between all the patches of size
    * (receptiveField x receptiveField)present in the image of the first layer.
    * This version of the method is not perfroming duplicated calculations.
    * Distance is pairwise operation and symetric. Hence after performing
    * operation for (x_1,x_2) there is no need for (x_2,x_1).
    * Only at the first layer calculation are performed on the on actual data.
    * @param m input matrix (image of the selected layer)
    * @param data storage for calculated distance measure for all the parts (output)
    * @param recFieldSize size of the receptive field for the current layer
    * @return 0 for success return
    */
    int calcDistanceMeasure(const Mat &m, vector<vector<double> > &data, int recFieldSize);

    /**
    * @brief Calculates distance between patches at the second layer and above
    *
    * Use this method together with precomputeDistance.
    * This function calculate difference between all the patches of size
    * (receptiveField x receptiveField) present in the image
    * at the current layer and returns it to vector of vectors.
    * This version of the method is perfroming duplicated calculations.
    * This is done for later convenience in segmentation function.
    * Distances are measured for the parts projected to the first layer.
    * @param m input matrix (image of the selected layer)
    * @param data storage for calculated distance measure for all the parts (output)
    * @param layerNum number of layer which is currently processed
    * @param recFieldSize size of the receptive field for the current layer
    * @return 0 for success return
    */
    int calcDistanceMeasureSec(const Mat &m,vector<vector<double> > &data, int layerNum, int recFieldSize);


    /**
    * @brief Precomputes the distance between parts on arbitrary layer.
    *
    * The method allows to obtain ceratain speed up in calculating the distance
    * between patches on the higher layers.
    * Distances are measured for the parts projected to the first layer.
    * Data are stored in the member variable layerDistance
    * @param layerNum number of layer which is currently processed
    * @param recFieldSize size of the receptive field for the current layer
    * @return 0 for success return
    */
    int precomputeDistance(int layerNum, int recFieldSize);

    /**
    * @brief Calculates receptive field with overlap.
    *
    * Receptive field at the botom layer is calculated as a median of values in
    * the neighbourhood specified by the receptive field size.
    * Border is replicated on the right and bottom part of the image.
    * @param m input matrix (image)
    * @param recept storage for an image after contraction
    * @param recFieldSize size of the receptive field
    * @return 0 for success return
    */
    int calcRecptFieldOverlap(const Mat &mOrient,const Mat &mMagnit, Mat &recept, int recFieldSize);

    /**
    * @brief Calculates normals using gradients calculated in X and Y direction.
    *
    * It is not the best way of calculating normals. Not whole information is
    * available since the projection from sphere into circle.
    * @param gradX input matrix with calculated values of gradient in X direction
    * @param gradY input matrix with calculated values of gradient in Y direction
    * @param normals storage for calculated normals
    * @return 0 for success return
    */
    int calcNormals(const Mat &gradX,const Mat &gradY, Mat &normals);
    /**
    * @brief Calculates normals using PCA.
    *
    * This method uses PCA to obtain the proper values of normals. The normals are
    * calculated using covariance matrix computed for the window on depth data.
    * The vector associated with the smallest eigen value is the normal.
    * @param m input matrix (image)
    * @param normals storage for calculated normals
    * @param recFieldSize size of the receptive field
    * @return 0 for success return
    */
    int calcNormalsPCA(const Mat &m, Mat &normals, int recFieldSize);

    /**
    * @brief Calculates distance measure between normals
    *
    * Use this method together with calcNormalsPCA().
    * This function calculate difference between all the patches of size
    * (receptiveField x receptiveField) present in the image and computed by
    * calcNormalsPCA(). Distance is obtained as a dot product between two patches.
    * The function is a resource consuming. While calculate all the possible
    * combintations between parts.
    * @param m input matrix (normals image)
    * @param data storage for calculated distance measure for all the parts (output)
    * @return 0 for success return
    */
    int calcDistanceMeasureNormals(const Mat &m, vector<vector<float> > &data);

    /**
    * @brief Calculates discrete versions of derivatives using Sobel mask.
    *
    * The method allows to obtain discrete values of the derivates of the input
    * image. Sobel mask was used, hence the scaling is performed after first derivative
    * is obtained and passed to the second derivative.
    * @param m input matrix (image)
    * @param gradX storage for calculated first derivative along X direction
    * @param gradY storage for calculated first derivative along Y direction
    * @param gradXX storage for calculated second derivative along X direction
    * @param gradYY storage for calculated second derivative along Y direction
    * @param gradXY storage for calculated mixed derivative first (X) second (Y)
    * @param kernelSize size of the kernel used for derivative
    * @return 0 for success return
    */
    int calcDerivatives(const Mat &m, Mat &gradX, Mat &gradY, Mat &gradXX, Mat &gradYY, Mat &gradXY,char kernelSize);
    /**
    * @brief Calculates discrete versions of local surface curvature descriptor.
    *
    * Use this method together with calcDerivatives().
    * The method allows to obtain Gaussian and mean curvature of the surfaces.
    * These are the local measures of surface shape. Due to the computation of
    * second derivativea they are very sensitive to noise. Intensive smoothing
    * is required.
    * @param m input matrix (image)
    * @param gradX input first derivative along X direction
    * @param gradY input first derivative along Y direction
    * @param gradXX input second derivative along X direction
    * @param gradYY input second derivative along Y direction
    * @param gradXY input mixed derivative first (X) second (Y)
    * @param H storage for computed values of mean curvature
    * @param K storage for computed values of Gaussian curvature
    * @param HS storage for computed signs of mean curvature
    * @param KS storage for computed signs of Gaussian curvature
    * @param HK storage for computed HK-map
    * @param kernelSize size of the kernel used for derivative
    * @return 0 for success return
    */
    int calcHKmap(const Mat &gradX, const Mat &gradY, const Mat &gradXX, const Mat &gradYY, const Mat &gradXY, Mat &H, Mat &K, Mat &HS, Mat &KS, Mat &HK);

    /**
    * @brief Loads a .ply file -- vertices and normals but not faces.
    *
    * Use this method together with calcDerivatives().
    * The method allows to obtain Gaussian and mean curvature of the surfaces.
    * These are the local measures of surface shape. Due to the computation of
    * second derivativea they are very sensitive to noise. Intensive smoothing
    * is required.
    * @param m input matrix (image)
    * @param gradX input first derivative along X direction
    * @param gradY input first derivative along Y direction
    * @param gradXX input second derivative along X direction
    * @param gradYY input second derivative along Y direction
    * @param gradXY input mixed derivative first (X) second (Y)
    * @param H storage for computed values of mean curvature
    * @param K storage for computed values of Gaussian curvature
    * @param HS storage for computed signs of mean curvature
    * @param KS storage for computed signs of Gaussian curvature
    * @param HK storage for computed HK-map
    * @param kernelSize size of the kernel used for derivative
    * @return 0 for success return
    */
    int readPlyFile(std::string fileName, std::vector<Eigen::Vector3f> &outputCloud, std::vector<Eigen::Vector3f> &outputNormals, std::vector<Eigen::Vector4i> &outputFaces);
    int readDiscretization(std::string fileName);
    int quantizeNormals(const std::vector<Eigen::Vector3f> &inputNormals, std::vector<Eigen::Vector3f> &outputNormals);
    int writePlyToFile(std::string fileName, const std::vector<Eigen::Vector3f> &inputCloud, const std::vector<Eigen::Vector3f> &inputNormals, const std::vector<Eigen::Vector4i> &inputFaces );

    void split(const string& s, char c,vector<string>& v);

    int obtainSecondLayerPartsHash(const Mat &m, int recFieldSize);
    int obtainSecondLayerParts(const Mat &m, int recFieldSize);
    int rarePartsRemoval(vector<vector<double> > &data, vector<vector<double> > &dataT, int layerNum);

    int computePartsFrequency(vector<vector<double> > &data, int layerNum, int difference);
    int reconstructImage(const Mat &m, vector<vector<double> > &data, Mat &output, int layerNum,int recFieldSize);
    int readROIs(path p,vector<Rect> &output);
    int readMultiImage(path p, vector<Mat> &output, vector<Rect> &input);

    ///Mapping parts with respective hash value
    multimap<long,Mat> _hashedParts;
    ///Mapping from Layer to Layer
    vector<unordered_map<int,Mat> >_layerMap;
    ///Frequency of each part at each layer
    vector<unordered_map<int,int> >_layerMapFrequency;

    ///Member variable for storring the precalculated distances for each layer
    ///_layerDistance[layerNumber][firstPartNumber][secondPartNumber] = distance
    vector<unordered_map<int, unordered_map<int,int> > > _layerDistance;

    ///Mapping a distance between normals in discretization function
    map<float,Eigen::Vector3f> _normalsDiscreteMap;



protected:

private:

    /**
    * @brief Calculate the median of the patch
    *
    * @param m input matrix (image)
    * @return position of the median pixel
    */
    int median(const Mat &src);

    /**
    * @brief Calculate the norm of two patches taking into account the circular
    * data.
    *
    * For example in angular degrees 10 deg and 350 deg are 20 deg away not 340.
    * @param a first input matrix (patch)
    * @param b second input matrix (patch)
    * @param N the lenght of the circular buffer
    * @param recFieldSize size of the receptive field
    * @return distance between two patches
    */
    float normCircular(const Mat &a,const Mat &b, int N, int recFieldSize);
    //Variables regarding gradient
    int ddepth_;
    int scale_;
    int delta_;
    int gaussSize_;

    //Variables regarding histogram
    ///Number of bins
    int histSize_;
    ///Range of values
    float range_[2];
    ///true -- same size of the bins
    bool uniform_;
    ///true -- clearing the histogram in the beginning
    bool accumulate_;
    ///Size of the image to display the histogram
    int histW_;
    int histH_;
    ///Number of Layers
    int numLayers;
    ///Size of the image for Layer_n
    vector<int> divImageXLayer;
    vector<int> divImageYLayer;
    ///Storage for parts on Layer_n
    vector<vector<Mat> > _parts;
    ///Storage for parts on Layer_n when using normals Vec3f
    vector<vector<Point3f> > _parts3f;

    ///Storing discretized version of Gaussian Sphere
    std::vector<Eigen::Vector3f> discretizedNormals;

};
