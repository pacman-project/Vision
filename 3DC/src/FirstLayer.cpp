#include "FirstLayer.h"



FirstLayer::FirstLayer(): ddepth_(CV_32F), scale_(1), delta_(0), gaussSize_(5),
    histSize_(360),uniform_(true),accumulate_(false), histW_(360),histH_(200),
    numLayers(5)
{
    range_[0] = -180;
    range_[1] = 180;
    divImageXLayer.resize(numLayers);
    divImageYLayer.resize(numLayers);
    _parts.resize(numLayers);
    _parts3f.resize(numLayers);
    _layerMap.resize(numLayers);
    _layerMapFrequency.resize(numLayers);
    _layerDistance.resize(numLayers);

}

int FirstLayer::matrixToTxtFile(const Mat &m, string filename){
    ofstream myFile;
    myFile.open (filename);
    myFile << "M (default) = " << endl <<     m       << endl << endl;
    myFile.close();

    return 0;
}

int FirstLayer::vectorToTxtFile(const vector<vector<double> > data, string filename){
    ofstream myFile;

    myFile.open (filename);
    myFile << "GLfloat vector[]={" << endl;
    for(int i = 0; i<data.size();i++){
        vector<double> row = data[i];
        copy(row.begin(), row.end(), ostream_iterator<double>(myFile , ", "));
        myFile << "}" << endl;
    }
   myFile.close();

    return 0;
}

int FirstLayer::vectorToTxtFile(const vector<vector<float> > data, string filename){
    ofstream myFile;

    myFile.open (filename);
    myFile << "GLfloat vector[]={" << endl;
    for(int i = 0; i<data.size();i++){
        vector<float> row = data[i];
        copy(row.begin(), row.end(), ostream_iterator<float>(myFile , ", "));
        myFile << "}" << endl;
    }
   myFile.close();

    return 0;
}

int FirstLayer::readMultiImage(path p,vector<Mat> &output, vector<Rect> &input){
    try
      {
        if (exists(p))    // does p actually exist?
        {
          if (is_regular_file(p))        // is p a regular file?
            cout << p << " size is " << file_size(p) << '\n';

          else if (is_directory(p))      // is p a directory?
          {
            cout << p << " is a directory containing:\n";

            typedef vector<path> vec;             // store paths,
            vec v;                                // so we can sort them later

            copy(directory_iterator(p), directory_iterator(), back_inserter(v));

            sort(v.begin(), v.end());             // sort, since directory iteration
                                                  // is not ordered on some file systems
            int counter = 0;
            for (vec::const_iterator it(v.begin()), it_end(v.end()); it != it_end; ++it)
            {
              cout << "   " << *it << '\n';
              Mat image = imread(it->string(), CV_LOAD_IMAGE_UNCHANGED);
              Mat imageCrop = image(input[counter]);
              output.push_back(imageCrop);
              counter++;
            }
          }
          else
            cout << p << " exists, but is neither a regular file nor a directory\n";
        }
        else
          cout << p << " does not exist\n";
      }

      catch (const filesystem_error& ex)
      {
        cout << ex.what() << '\n';
      }
}

int FirstLayer::readROIs(path p,vector<Rect> &output){
    try
      {
        if (exists(p))    // does p actually exist?
        {
          if (is_regular_file(p))        // is p a regular file?
            cout << p << " size is " << file_size(p) << '\n';

          else if (is_directory(p))      // is p a directory?
          {
            cout << p << " is a directory containing:\n";

            typedef vector<path> vec;             // store paths,
            vec v;                                // so we can sort them later

            copy(directory_iterator(p), directory_iterator(), back_inserter(v));

            sort(v.begin(), v.end());             // sort, since directory iteration
                                                  // is not ordered on some file systems

            for (vec::const_iterator it(v.begin()), it_end(v.end()); it != it_end; ++it)
            {
              cout << "   " << *it << '\n';
              string line;
                ifstream myfile (it->string());
                if (myfile.is_open())
                {
                  while ( getline (myfile,line) )
                  {
                    stringstream linestream(line);
                    vector<string> strNumber;
                    for(int i = 0; i < 5; i++){
                        string tempString;
                        getline(linestream, tempString, ',');
                        strNumber.push_back(tempString);
                    }
                    output.push_back(Rect(stoi(strNumber[0]),stoi(strNumber[1]),stoi(strNumber[2]),stoi(strNumber[3])));

                  }
                  myfile.close();
                }

                else cout << "Unable to open file";
            }
          }
          else
            cout << p << " exists, but is neither a regular file nor a directory\n";
        }
        else
          cout << p << " does not exist\n";
      }

      catch (const filesystem_error& ex)
      {
        cout << ex.what() << '\n';
      }
}

int FirstLayer::calcGradient(const Mat &m, Mat &gradX, Mat &gradY, Mat &gradOrient, Mat &gradMagnit){
    Mat img;
    //m.copyTo(img);
    // Smoothing
    GaussianBlur( m, img, Size(gaussSize_,gaussSize_), 0, 0, BORDER_DEFAULT );
    // Gradient X
    Sobel(img, gradX, ddepth_, 1, 0, CV_SCHARR, scale_, delta_, BORDER_DEFAULT );
    // Gradient Y
    Sobel(img, gradY, ddepth_, 0, 1, CV_SCHARR, scale_, delta_, BORDER_DEFAULT);
    //Calculate direction of gradient
    gradOrient.create(gradX.rows,gradX.cols,ddepth_);
    phase(gradX,gradY,gradOrient,true);
    magnitude(gradX,gradY,gradMagnit);
    //gradient was scalled
    //gradMagnit = gradMagnit*0.01;
    //cout << gradMagnit << endl ;

//    for(int i=0; i<gradX.cols; i++){
//              for(int j=0; j<gradX.rows; j++){
//                 gradOrient.at<double>(j,i) = 180/3.14*(atan2( gradY.at<double>(j,i),gradX.at<double>(j,i)));
//              }
//    }

    return 0;
}


int FirstLayer::calcDerivatives(const Mat &m, Mat &gradX, Mat &gradY, Mat &gradXX, Mat &gradYY, Mat &gradXY, char kernelSize){
    Mat img;
    //m.copyTo(img);
    // Smoothing
    GaussianBlur( m, img, Size(kernelSize,kernelSize), 3*kernelSize, 3*kernelSize, BORDER_DEFAULT );
    //bilateralFilter(m,img,kernelSize,kernelSize*2,kernelSize/2);
    //scalling of first order derivatives according to kernel size
    double scaleFactor =  pow(double(kernelSize),double(4))*pow(double(10),double((-kernelSize)));

    // Gradient X
    Sobel(img, gradX, ddepth_, 1, 0, kernelSize, scaleFactor, delta_, BORDER_DEFAULT );
    // Gradient Y
    Sobel(img, gradY, ddepth_, 0, 1, kernelSize, scaleFactor, delta_, BORDER_DEFAULT);
    //Gradient XX
    Sobel(gradX, gradXX, ddepth_, 1, 0, kernelSize, scale_, delta_, BORDER_DEFAULT );
    // Gradient YY
    Sobel(gradY, gradYY, ddepth_, 0, 1, kernelSize, scale_, delta_, BORDER_DEFAULT);
    // Gradient XY
    Sobel(gradX, gradXY, ddepth_, 0, 1, kernelSize, scale_, delta_, BORDER_DEFAULT );



    return 0;
}

int FirstLayer::calcHKmap(const Mat &gradX, const Mat &gradY, const Mat &gradXX, const Mat &gradYY, const Mat &gradXY, Mat &H, Mat &K, Mat &HS, Mat &KS, Mat &HK){
    cout << "Size of gradX: " << gradX.cols <<" x " << gradX.rows << endl;
    H.create(gradX.rows,gradX.cols,ddepth_);
    K.create(gradX.rows,gradX.cols,ddepth_);
    HK.create(gradX.rows,gradX.cols,CV_8UC1);
    HS.create(gradX.rows,gradX.cols,CV_8SC1);
    KS.create(gradX.rows,gradX.cols,CV_8SC1);

    for(int i=0; i<gradX.rows; i++){
          for(int j=0; j<gradX.cols; j++){
              float gXX = gradXX.at<float>(i,j);
              float gYY = gradYY.at<float>(i,j);
              float gXY = gradXY.at<float>(i,j);
              float gX = gradX.at<float>(i,j);
              float gY = gradY.at<float>(i,j);
              float valK = (gXX*gYY-(gXY*gXY))/((1+gX*gX+gY*gY)*(1+gX*gX+gY*gY));
              float valH = -0.5*((1+gY*gY)*gXX+(1+gX*gX)*gYY-2*gX*gY*gXY)/
                      sqrt((1+gX*gX+gY*gY)*(1+gX*gX+gY*gY)*(1+gX*gX+gY*gY));
              K.at<float>(i,j) = valK;
              H.at<float>(i,j) = valH;
              char signK = 0;
              char signH = 0;
              //tolerance parameter
              int epsilon = 50;
              //related to size of the derivative kernel
              //related to the variance in second derivative or in K or H
              if (valK > epsilon) signK = 1;
              if (valK < -epsilon) signK = -1;
              if (valH > epsilon) signH = 1;
              if (valH < -epsilon) signH = -1;
              HS.at<char>(i,j) = signH;
              KS.at<char>(i,j) = signK;
              HK.at<char>(i,j) = 1 + 3*(1+signH)+(1-signK);
          }
      }

    //GaussianBlur( HKtemp, HK, Size(gaussSize_,gaussSize_), 5, 5, BORDER_DEFAULT );
    return 0;
}

int FirstLayer::calcNormals(const Mat &gradX,const Mat &gradY, Mat &normals){

    //Calculate direction of gradient
    normals.create(gradX.rows,gradX.cols,CV_32FC3);


    for(int i=0; i<gradX.cols; i++){
              for(int j=0; j<gradX.rows; j++){
                 float normX = -gradX.at<float>(j,i);
                 float normY = -gradY.at<float>(j,i);
                 float normZ = 1;
                 float len = sqrt(normX*normX+normY*normY+1);
                 normals.at<cv::Vec3f>(j,i)[0]= normX/len;
                 normals.at<cv::Vec3f>(j,i)[1]= normY/len;
                 normals.at<cv::Vec3f>(j,i)[2]= normZ/len;
                 }
    }
    return 0;
}


//This function in equivalent to calc recept field
//Look for better normals estimation
//Need to do this with overlap
int FirstLayer::calcNormalsPCA(const Mat &m, Mat &normals, int recFieldSize){

    Mat mT;
    copyMakeBorder(m,mT,0,(recFieldSize-m.rows%recFieldSize),0,(recFieldSize-m.cols%recFieldSize),BORDER_REPLICATE);
    divImageXLayer[0] = mT.cols/recFieldSize;
    divImageYLayer[0] = mT.rows/recFieldSize;
    //cout << "divs of image layer 0:" << divImageXLayer[0] << divImageYLayer[0];
    //cout << "A=[";
    normals.create(divImageXLayer[0],divImageYLayer[0],CV_32FC3);
    for(int i=0; i<divImageXLayer[0]; i++){
          for(int j=0; j<divImageYLayer[0]; j++){
              Mat subImg = mT(Rect(i*recFieldSize, j*recFieldSize,recFieldSize, recFieldSize));
              PCA pca(subImg, // pass the data
                         Mat(), // there is no pre-computed mean vector,
                                // so let the PCA engine to compute it
                         CV_PCA_DATA_AS_ROW, // indicate that the vectors
                                             // are stored as matrix rows
                                             // (use CV_PCA_DATA_AS_COL if the vectors are
                                             // the matrix columns)
                         6 // specify how many principal components to retain
                         );
              normals.at<cv::Vec3f>(i,j) = pca.eigenvectors.row(2);
              //cout << pca.eigenvectors.row(2) <<";"<< endl;
              //cout << pca.eigenvalues << endl;
          }
      }
    return 0;
}


int FirstLayer::calcDistanceMeasureNormals(const Mat &m, vector<vector<float> > &data){
    for(int i=0; i<1; i++){
          for(int j=0; j<1; j++){
              Point3f subImg = m.at<cv::Vec3f>(i,j);
              _parts3f[0].push_back(subImg);
              vector<float> row;
              for(int k=0; k<divImageXLayer[0]; k++){
                    for(int l=0; l<divImageYLayer[0]; l++){
                        Point3f subImgTemp = m.at<cv::Vec3f>(k,l);
                        row.push_back(acos(subImg.dot(subImgTemp))*180/3.14159265);
                    }
              }
              data.push_back(row);
          }
      }
    return 0;

}

int FirstLayer::calcHistogram(const Mat &m, Mat &hist){
    Mat img;
    m.convertTo(img,CV_32F);
    const float* histRange = { range_ };
    Mat depthHist;
    calcHist( &img, 1, 0, Mat(), depthHist, 1, &histSize_, &histRange, uniform_, accumulate_ );
    int bin_w = cvRound( (double) histW_/histSize_ );
    hist.create(histH_, histW_, CV_8UC3);
    normalize(depthHist, depthHist, 0, hist.rows, NORM_MINMAX, -1, Mat() );
    for( int i = 1; i < histSize_; i++ )
    {
      line( hist, Point( bin_w*(i-1), histH_ - cvRound(depthHist.at<float>(i-1)) ) ,
                       Point( bin_w*(i), histH_ - cvRound(depthHist.at<float>(i)) ),
                       Scalar( 255, 0, 0), 2, 8, 0  );
    }

    return 0;
}

int FirstLayer::calcColorMap(const Mat &m, Mat &colorMap){

    double min;
    double max;
    minMaxIdx(m, &min, &max);
    cout << "It is min, max: " << min << "; " << max << endl;
    cv::Mat adjMap;
    cv::convertScaleAbs(m, adjMap, 255 / max);
    //applyColorMap(adjMap, colorMap, cv::COLORMAP_PINK);  // good for gradients
    applyColorMap(adjMap, colorMap, cv::COLORMAP_HSV);  // good for orientation
    return 0;
}

//assuming odd number of elements and square shape receptive field
int FirstLayer::median(const Mat &src){
    vector<float> data;
    int iter = (src.cols+1)/2;
    for(int i=0; i<src.cols; i++){
        for(int j=0; j<src.rows; j++){
            data.push_back(src.at<float>(i,j));
        }
    }
    return data[iter];
}

int FirstLayer::calcRecptFieldOverlap(const Mat &mOrient,const Mat &mMagnit, Mat &recept, int recFieldSize){
    Mat mOrientB;
    medianBlur(mOrient,mOrientB,3);  //after blur
    Mat mOrientT;

    copyMakeBorder(mOrientB,mOrientT,1,(recFieldSize-mOrient.rows%recFieldSize)+1,1,(recFieldSize-mOrient.cols%recFieldSize)+1,BORDER_REPLICATE);
    Mat mMagnitT;
    copyMakeBorder(mMagnit,mMagnitT,1,(recFieldSize-mMagnit.rows%recFieldSize)+1,1,(recFieldSize-mMagnit.cols%recFieldSize)+1,BORDER_REPLICATE);

    divImageXLayer[0] = mOrientT.rows/recFieldSize;
    divImageYLayer[0] = mOrientT.cols/recFieldSize;

    recept.create(divImageXLayer[0],divImageYLayer[0],CV_32S);
    vector<int> parts;
    for(int i=0; i<divImageXLayer[0]; i++){
          for(int j=0; j<divImageYLayer[0]; j++){
//aimed at using 3X3 recField size, larger recField requrires testing
              vector<double> stdMagnitValues;
              vector<double> orientValues;
              for(int ii=-1; ii<recFieldSize-1; ii++){
                  for(int jj=-1; jj<recFieldSize-1; jj++){
                      Scalar     mean;
                      Scalar     stddev;
                      Mat subImgO = mOrientT(Rect((j*recFieldSize)+jj+1, (i*recFieldSize)+ii+1,recFieldSize, recFieldSize));
                      orientValues.push_back(subImgO.at<float>((recFieldSize+1)/2,(recFieldSize+1)/2));
                      Mat subImgM = mMagnitT(Rect((j*recFieldSize)+jj+1, (i*recFieldSize)+ii+1,recFieldSize, recFieldSize));
                      cv::meanStdDev (subImgM, mean, stddev );
                      stdMagnitValues.push_back(stddev.val[0]);
                  }
              }
              int minPos = distance(stdMagnitValues.begin(),min_element(stdMagnitValues.begin(),stdMagnitValues.end()));
              int medValue = orientValues[minPos];
              recept.at<int>(i,j) = medValue;
              parts.push_back(medValue);
          }
      }
    std::sort (parts.begin(), parts.end());
    std::vector<int>::iterator it;
    it = std::unique (parts.begin(), parts.end());
    parts.resize( std::distance(parts.begin(),it) );
    std::cout << "Number of cols and rows at layer 1: " << divImageXLayer[0] << "; " << divImageYLayer[0] << std::endl;
    cout << "Number of parts at layer 1ab: " << parts.size() << endl;
    return 0;
}


int FirstLayer::calcRecptField(const Mat &m, Mat &recept, int recFieldSize){

    Mat mT;
    copyMakeBorder(m,mT,0,(recFieldSize-m.rows%recFieldSize),0,(recFieldSize-m.cols%recFieldSize),BORDER_REPLICATE);
    divImageXLayer[0] = mT.rows/recFieldSize;
    divImageYLayer[0] = mT.cols/recFieldSize;

    recept.create(divImageXLayer[0],divImageYLayer[0],CV_32S);
    vector<int> parts;
    for(int i=0; i<divImageXLayer[0]; i++){
          for(int j=0; j<divImageYLayer[0]; j++){
              Mat subImg = mT(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
//using mean
//              cv:Scalar tempVal = mean( subImg );
//              recept.at<int>(i,j) = tempVal.val[0];
//              parts.push_back(tempVal.val[0]);
//using median
              int medValue = median(subImg);
              recept.at<int>(i,j) = medValue;
              parts.push_back(medValue);
          }
      }
    std::sort (parts.begin(), parts.end());
    std::vector<int>::iterator it;
    it = std::unique (parts.begin(), parts.end());
    parts.resize( std::distance(parts.begin(),it) );
    std::cout << "Number of cols and rows at layer 1: " << divImageXLayer[0] << "; " << divImageYLayer[0] << std::endl;
    cout << "Number of parts at layer 1: " << parts.size() << endl;
    return 0;
}

int FirstLayer::obtainSecondLayerParts(const Mat &m, int recFieldSize){
    Mat mT;
    copyMakeBorder(m,mT,1,(recFieldSize-m.rows%recFieldSize)+1,1,(recFieldSize-m.cols%recFieldSize)+1,BORDER_REPLICATE);

    divImageXLayer[1] = mT.rows/recFieldSize;
    divImageYLayer[1] = mT.cols/recFieldSize;

    std::cout << "Number of cols and rows at layer 2: " << divImageXLayer[1] << "; " << divImageYLayer[1] << std::endl;

    for(int i=0; i<divImageXLayer[1]; i++){
          for(int j=0; j<divImageYLayer[1]; j++){
              Mat subImg = mT(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
              //_parts[1].push_back(subImg); //will be stored in multimap
              Scalar     mean;
              Scalar     stddev;
              bitset<16> differences;
              int       meanVal;
              cv::meanStdDev (subImg, mean, stddev);
              meanVal = int(mean.val[0]);
              //computing averaged hash function
              int counter = 0;
              for(int ii=0; ii<recFieldSize; ii++){
                  for(int jj=0; jj<recFieldSize; jj++){
                      if( subImg.at<int>(ii,jj) > meanVal) differences[counter] = 1;
                      else differences[counter] = 0;
                      counter++;
                  }
              }
              _hashedParts.insert(std::pair<long,Mat>(differences.to_ulong(),subImg));
           }
    }
    //here starts processing of the hashed map
    vector<pair<long,Mat> > myVector;
    vector<long> myVectorLong;
    //each key(bucket) can contain several patches, have to check if they are the same
    //what's the resolution of hashing?
    //this trick is used to obtain the unique keys in the multimap
    unique_copy(begin(_hashedParts),
                  end(_hashedParts),
                  back_inserter(myVector),
                  [](const pair<long,Mat> &entry1,
                     const pair<long,Mat> &entry2) {
                       return (entry1.first == entry2.first);
                   }
                 );
    //**********************
    //obtaining hash values only
    for(auto it = myVector.begin(); it != myVector.end(); ++it)
        myVectorLong.push_back(it->first);
    cout << "Number of parts in hashMultiMap: " << _hashedParts.size() << endl;
    cout << "Number of parts in hashMultiMapVector: " << myVectorLong.size() << endl;
    int countRelevant = 0;
    for (auto it = myVectorLong.begin(); it != myVectorLong.end(); it++){
           int numParts = _hashedParts.count(*it);
           cout << "Number of elements in bucket " << *it << ": "<<  numParts << endl;
               multimap<long,Mat>::iterator pos;
               cout << *it << ": " << endl;
               for (pos = _hashedParts.lower_bound(*it);
                    pos != _hashedParts.upper_bound(*it); ++pos) {
                       cout << "    " << pos->second << endl;
               }
           if (numParts > 1){
                  countRelevant++;
           }
    }
    cout << "Number of relevant parts: " << countRelevant << endl;
    //sort(myVectorLong.begin(),myVectorLong.end()); //it is not really necesary
    //building a histogram
//    int bucket_size = 10;
//    int number_of_buckets = (int)ceil(*max_element(myVectorLong.begin(),myVectorLong.end()) / bucket_size);
//    cout << "Number of buckets: " << number_of_buckets << endl;
//    vector<int> histogram(number_of_buckets);
//    for(auto it = myVectorLong.begin(); it != myVectorLong.end(); ++it){
//            int bucket = (int)floor(*it / bucket_size);
//            //if something fails perform range check
//            histogram[bucket] += 1;
//    }
//    int counter = 0;
//    vector<int> availableHashes;
//    for(auto it = histogram.begin(); it != histogram.end(); ++it){
//        if(*it > 1){  //resolve the issue with number of occurences
//            availableHashes.push_back(counter+5);
//            cout << *it << "; ";
//        }
//        counter++;
//    }
//    cout << endl;
//    cout << "Number of elements at the end: " << availableHashes.size() << endl;

    return 0;

}

//int FirstLayer::obtainSecondLayerParts(const Mat &m, int recFieldSize){
//    Mat mT;
//    copyMakeBorder(m,mT,1,(recFieldSize-m.rows%recFieldSize)+1,1,(recFieldSize-m.cols%recFieldSize)+1,BORDER_REPLICATE);

//    divImageXLayer[1] = mT.rows/recFieldSize;
//    divImageYLayer[1] = mT.cols/recFieldSize;

//    std::cout << "Number of cols and rows at layer 2: " << divImageXLayer[1] << "; " << divImageYLayer[1] << std::endl;

//    for(int i=0; i<divImageXLayer[1]; i++){
//          for(int j=0; j<divImageYLayer[1]; j++){
//              Mat subImg = mT(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
//              _parts[1].push_back(subImg);
//           }
//    }
//    return 0;
//}



int FirstLayer::calcDistanceMeasure(const Mat &m, vector<vector<double> > &data, int recFieldSize){
    divImageXLayer[1] = m.rows/recFieldSize;
    divImageYLayer[1] = m.cols/recFieldSize;
    std::cout << "Number of cols and rows at layer 2: " << divImageXLayer[1] << "; " << divImageYLayer[1] << std::endl;

    for(int i=0; i<divImageXLayer[1]; i++){
          for(int j=0; j<divImageYLayer[1]; j++){
              Mat subImg = m(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
              _parts[1].push_back(subImg);
              vector<double> row;
              for(int k=0; k<divImageXLayer[1]; k++){
                    for(int l=0; l<divImageYLayer[1]; l++){
                        Mat subImgTemp = m(Rect(l*recFieldSize, k*recFieldSize,recFieldSize, recFieldSize));
                        row.push_back(norm(subImg,subImgTemp));
                    }
                  }
              data.push_back(row);
          }
      }
    cout << "Size of parts: " << _parts[1].size() << endl;
    return 0;

}


int FirstLayer::calcDistanceMeasureHalf(const Mat &m, vector<vector<double> > &data, int recFieldSize){
    divImageXLayer[1] = m.rows/recFieldSize;
    divImageYLayer[1] = m.cols/recFieldSize;
    std::cout << "Number of cols and rows at layer 2: " << divImageXLayer[1] << "; " << divImageYLayer[1] << std::endl;

    for(int i=0; i<divImageXLayer[1]; i++){
          for(int j=0; j<divImageYLayer[1]; j++){
              Mat subImg = m(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
              _parts[1].push_back(subImg);
              vector<double> row;
              for(int k=i; k<divImageXLayer[1]; k++){
                  if(k == i){
                      for(int l=j; l<divImageYLayer[1]; l++){
                          Mat subImgTemp = m(Rect(l*recFieldSize, k*recFieldSize,recFieldSize, recFieldSize));
                          row.push_back(norm(subImg,subImgTemp));
                      }
                  }
                  else{
                    for(int l=0; l<divImageYLayer[1]; l++){
                        Mat subImgTemp = m(Rect(l*recFieldSize, k*recFieldSize,recFieldSize, recFieldSize));
                        row.push_back(norm(subImg,subImgTemp));
                    }
                  }
              }
              data.push_back(row);
          }
      }
    return 0;

}

int FirstLayer::calcDistanceMeasureSec(const Mat &m,vector<vector<double> > &data,  int layerNum, int recFieldSize){
    divImageXLayer[layerNum-1] = m.rows/recFieldSize;
    divImageYLayer[layerNum-1] = m.cols/recFieldSize;
    std::cout << "Number of cols and rows at layer " << layerNum<<": " << divImageXLayer[layerNum-1] << "; " << divImageYLayer[layerNum-1] << std::endl;

    for(int i=0; i<divImageXLayer[layerNum-1]; i++){
          for(int j=0; j<divImageYLayer[layerNum-1]; j++){
              Mat subImg = m(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
              _parts[layerNum-1].push_back(subImg);
              vector<double> row;
              for(int k=0; k<divImageXLayer[layerNum-1]; k++){
                    for(int l=0; l<divImageYLayer[layerNum-1]; l++){
                        Mat subImgTemp = m(Rect(l*recFieldSize, k*recFieldSize,recFieldSize, recFieldSize));
                        int cumulateNorm = 0;
                        for(int ii=0; ii<recFieldSize; ii++){
                            for(int jj=0; jj<recFieldSize; jj++){
                                int a = subImg.at<int>(ii,jj);
                                int b = subImgTemp.at<int>(ii,jj);
                                //cumulateNorm+=(norm(_layerMap[1][a],_layerMap[1][b]));
                                cumulateNorm+=_layerDistance[layerNum-2][a][b];
                            }
                        }
                        row.push_back(cumulateNorm);
                    }
              }
              data.push_back(row);
          }
      }
    return 0;
}


int FirstLayer::segmentParts(vector<vector<double> > &data, Mat &output, int layerNum, int difference){
    output.create(divImageXLayer[layerNum-1],divImageYLayer[layerNum-1],CV_32S);
    for(int i = 0; i<data.size();i++){
        vector<double> row = data[i];
        int partsCount = 0;
        if(!(std::isnan(row[0]))){
            for(int j=0; j<row.size();j++){
                if(row[j]<difference){
                    output.at<int>(j/divImageYLayer[layerNum-1],j%divImageYLayer[layerNum-1]) = i;
                    partsCount++;
                    data[j][0]=std::numeric_limits<double>::quiet_NaN();
                }
            }
            //moved from inside the row[j]<difference
        _layerMap[layerNum-1][i] = _parts[layerNum-1][i]; //writing here several times
        _layerMapFrequency[layerNum-1][i] = partsCount; //writing here several times
        }
    }
    return 0;

//      std::cout << "mymap's buckets contain:\n";
//        for ( unsigned ii = 0; ii < layerMap[layerNum-1].bucket_count(); ++ii) {
//          std::cout << "bucket #" << ii << " contains:";
//          for ( auto local_it = layerMap[layerNum-1].begin(ii); local_it!= layerMap[layerNum-1].end(ii); ++local_it )
//            std::cout << " " << local_it->first << ":" << local_it->second;
//          std::cout << std::endl;
//        }
}

int FirstLayer::reconstructImage(const Mat &m, vector<vector<double> > &data, Mat &output, int layerNum,int recFieldSize){
    output.create(divImageXLayer[layerNum-1],divImageYLayer[layerNum-1],CV_32S);
    for(int i=0; i<divImageXLayer[layerNum-1]; i++){
          for(int j=0; j<divImageYLayer[layerNum-1]; j++){
              Mat subImg = m(Rect(j*recFieldSize, i*recFieldSize,recFieldSize, recFieldSize));
              vector<double> distances;
              vector<int> relatedParts;
              for ( auto it = _layerMap[layerNum-1].begin(); it != _layerMap[layerNum-1].end(); ++it ){
                  distances.push_back(normCircular(it->second,subImg,360,recFieldSize));
                  relatedParts.push_back(it->first);
              }
              int minPos = distance(distances.begin(),min_element(distances.begin(),distances.end()));
              output.at<int>(i,j) = relatedParts.at(minPos);
           }
      }
    return 0;
}

int FirstLayer::computePartsFrequency(vector<vector<double> > &data, int layerNum, int difference){
    for(int i = 0; i<data.size();i++){
        vector<double> row = data[i];
        int partsCount = 0;
        if(!(std::isnan(row[0]))){
            for(int j=0; j<row.size();j++){
                if(row[j]<difference){
                    partsCount++;
                    if(i != j)
                    data[j][0]=std::numeric_limits<double>::quiet_NaN();
                }
            }
        _layerMap[layerNum-1][i] = _parts[layerNum-1][i]; //writing here several times
        _layerMapFrequency[layerNum-1][i] = partsCount; //writing here several times
        }
    }
    return 0;
}

int FirstLayer::rarePartsRemoval(vector<vector<double> > &dataT, vector<vector<double> > &data, int layerNum){
    std::cout << "mymap contains:";

    vector<int> rareParts;
        //acctually removing parts from layerMap
      for ( auto it = _layerMapFrequency[layerNum-1].begin(); it != _layerMapFrequency[layerNum-1].end(); ++it ){
        std::cout << " " << it->first << ":" << it->second;
        //how many occurrences to remove
        if (it->second <2){
            _layerMap[layerNum-1].erase(it->first);
            rareParts.push_back(it->first);
        }
      }
      cout << endl;
      //acctually removing parts from layerMapFrequency it couldn't be done in previous loop
      //while layerMapFrequency was used as an iterator
      for( auto it = rareParts.begin(); it != rareParts.end();it++)
          _layerMapFrequency[layerNum-1].erase(*it);

      //here parts will be exchanged
      for( auto it = rareParts.begin(); it != rareParts.end();it++){
        vector<double> row = dataT[*it];
        vector<double> distances;
        vector<double> relatedParts;
        for ( auto itt = _layerMapFrequency[layerNum-1].begin(); itt != _layerMapFrequency[layerNum-1].end(); ++itt ){
                distances.push_back(row[itt->first]);
                relatedParts.push_back(itt->first);
        }
        //cout << "Min element:" << *min_element(distances.begin(),distances.end());
        int minPos = distance(distances.begin(),min_element(distances.begin(),distances.end()));
        //cout << "  Part selected: " << relatedParts.at(minPos) << endl;
        //changing the distance matrix
        data[*it][0] = std::numeric_limits<double>::quiet_NaN();
        data[relatedParts.at(minPos)][*it] = 0;
        _layerMapFrequency[layerNum-1][relatedParts.at(minPos)]++;
      }

      std::cout << std::endl;
      std::cout << std::endl;
   std::cout << "mymap2 contains:";
   for ( auto it = _layerMapFrequency[layerNum-1].begin(); it != _layerMapFrequency[layerNum-1].end(); ++it )
        std::cout << " " << it->first << ":" << it->second;
    std::cout << std::endl;
    std::cout << std::endl;
}

//look at the datatypes
float FirstLayer::normCircular(const Mat &Ma,const Mat &Mb, int N, int recFieldSize){

    int norm = 0;
    for(int ii=0; ii<recFieldSize; ii++){
        for(int jj=0; jj<recFieldSize; jj++){
            int a = Ma.at<int>(ii,jj);
            int b = Mb.at<int>(ii,jj);
            int dm = int(abs(a - b))%N;
            if (dm < (N/2))
                norm += dm*dm;
            else
                norm += (N - dm)*(N - dm);

       }
    }
    return sqrt(norm);

}



//find an appropriate matching between the layers in precalculation and maps
int FirstLayer::precomputeDistance(int layerNum=2, int recFieldSize=3){
//think about circular distances  http://www.approxion.com/?p=337
    if(layerNum == 2){
        for ( auto it = _layerMap[layerNum-1].begin(); it != _layerMap[layerNum-1].end(); ++it ){
            for ( auto itt = _layerMap[layerNum-1].begin(); itt != _layerMap[layerNum-1].end(); ++itt ){
                //_layerDistance[layerNum-1][it->first][itt->first] = norm(it->second,itt->second);
                _layerDistance[layerNum-1][it->first][itt->first] = normCircular(it->second,itt->second,360, recFieldSize);
            }
        }
    }
    if(layerNum>2){
        for ( auto it = _layerMap[layerNum-1].begin(); it != _layerMap[layerNum-1].end(); ++it ){
            for ( auto itt = _layerMap[layerNum-1].begin(); itt != _layerMap[layerNum-1].end(); ++itt ){
                int cumulateNorm = 0;
                for(int ii=0; ii<recFieldSize; ii++){
                    for(int jj=0; jj<recFieldSize; jj++){
                        int a = it->second.at<int>(ii,jj);
                        int b = itt->second.at<int>(ii,jj);
                        cumulateNorm+=_layerDistance[layerNum-2][a][b];
                    }
                }
                _layerDistance[layerNum-1][it->first][itt->first] = cumulateNorm;
            }
        }

    }
    return 0;
}

void FirstLayer::split(const string& s, char c,vector<string>& v) {
   int i = 0;
   int j = s.find(c);

   while (j >= 0) {
      v.push_back(s.substr(i, j-i));
      i = ++j;
      j = s.find(c, j);

      if (j < 0) {
         v.push_back(s.substr(i, s.length()));
      }
   }
}

int FirstLayer::readPlyFile(std::string fileName, std::vector<Eigen::Vector3f> &outputCloud, std::vector<Eigen::Vector3f> &outputNormals, std::vector<Eigen::Vector4i> &outputFaces){
    std::cout << "Loading file: " << fileName << std::endl;
    std::string line;
    std::ifstream myfile (fileName.c_str());
    int j = 0;
    unsigned long vertexElements;
    unsigned long faceElements;
    if (myfile.is_open()){
        while ( getline (myfile,line) ){
            std::vector<string> lineSeq;
            split(line,' ',lineSeq);
           if(lineSeq.size() > 0){
                if(lineSeq[0] == "element"){
                    if(lineSeq[1] == "vertex"){
                        vertexElements = ::atof(lineSeq[2].c_str());
                        std::cout << "Number of vertices: " << vertexElements <<std::endl;
                    }
                    if(lineSeq[1] == "face"){
                        faceElements = ::atof(lineSeq[2].c_str());
                        std::cout << "Number of faces: " << faceElements <<std::endl;
                    }

                }
            }
            if(line == "end_header"){
                break;
            }
        }
        for(unsigned long  i=0; i < vertexElements; i++ ){
        getline (myfile,line);
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
            std::vector<string> lineSeq;
            split(line,' ',lineSeq);
            Eigen::Vector3f vecPoint;
            Eigen::Vector3f vecNorm;
            for(int i = 0; i< 3; i++) vecPoint(i) = ::atof(lineSeq[i].c_str());
            for(int i = 3; i< 6; i++) vecNorm(i-3) = ::atof(lineSeq[i].c_str());
            outputCloud.push_back(vecPoint);
            outputNormals.push_back(vecNorm);
        }
        std::cout << "Points and normals successfuly loaded from " << fileName << std::endl;
        std::cout << "The data occupies " << (2*(sizeof(std::vector<Eigen::Vector3f>) + (sizeof(Eigen::Vector3f) * outputCloud.size())))/1024/1024 << " MB" << std::endl;
        for(unsigned long  i=0; i < faceElements; i++ ){
        getline (myfile,line);
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
            std::vector<string> lineSeq;
            split(line,' ',lineSeq);
            Eigen::Vector4i vecFace;
            for(int i = 1; i< 5; i++) vecFace(i-1) = ::atof(lineSeq[i].c_str());
            outputFaces.push_back(vecFace);
        }
        myfile.close();

    }
    else{
        std::cout << "Unable to open ply file";
        return 1;
    }
    return 0;
}

int FirstLayer::readDiscretization(std::string fileName){
    std::cout << "Loading file: " << fileName << std::endl;
    std::string line;
    std::ifstream myfile (fileName.c_str());
    if (myfile.is_open()){
        while ( getline (myfile,line) ){
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
            std::vector<string> lineSeq;
            split(line,',',lineSeq);
            Eigen::Vector3f vecNorm;
            for(int i = 0; i< 3; i++) vecNorm(i) = ::atof(lineSeq[i].c_str());
            discretizedNormals.push_back(vecNorm);
        }
        std::cout << "Normals successfuly loaded from " << fileName << std::endl;
        std::cout << "The data occupies " << (sizeof(std::vector<Eigen::Vector3f>) + (sizeof(Eigen::Vector3f) * discretizedNormals.size()))/1024 << " kB" << std::endl;
        myfile.close();
    }
    else{
        std::cout << "Unable to open ply file";
        return 1;
    }
    for(int i = 0; i < discretizedNormals.size(); i++){
        float dist;
        float dot;
        dot = discretizedNormals[0].dot(discretizedNormals[i]);
        if(dot > 1) dot = 1;
        if(dot < -1) dot = -1;
        dist = acos(dot);
        cout << discretizedNormals[0].dot(discretizedNormals[i]) << std::endl;

        _normalsDiscreteMap[dist] = discretizedNormals[i];
    }
    std::cout << "mymap contains:";
      for ( auto it = _normalsDiscreteMap.begin(); it != _normalsDiscreteMap.end(); ++it )
        std::cout << " " << it->first << ":" << std::endl << it->second << std::endl;
        std::cout << std::endl;
    return 0;
}

int FirstLayer::quantizeNormals(const std::vector<Eigen::Vector3f> &inputNormals, std::vector<Eigen::Vector3f> &outputNormals){
    std::cout << "We are quantizing" << std::endl;
    for(unsigned long i =0; i < inputNormals.size(); i++ ){
        float dist;
        float dot;
        dot = discretizedNormals[0].dot(inputNormals[i]);
        if(dot > 1) dot = 1;
        if(dot < -1) dot = -1;
        dist = acos(dot);
    auto range = _normalsDiscreteMap.lower_bound(dist);
    outputNormals.push_back(range->second);
    }
    return 0;
}


int FirstLayer::writePlyToFile(std::string fileName, const std::vector<Eigen::Vector3f> &inputCloud, const std::vector<Eigen::Vector3f> &inputNormals, const std::vector<Eigen::Vector4i> &inputFaces ){
    ofstream myfile;
    myfile.open(fileName);
    myfile << "ply\r\n";
    myfile << "format ascii 1.0\r\n";
    myfile << "comment file created with Full3DHierearchy\r\n";
    myfile << "element vertex " << inputCloud.size() << "\r\n";
    myfile << "property float x\r\n";
    myfile << "property float y\r\n";
    myfile << "property float z\r\n";
    myfile << "property float nx\r\n";
    myfile << "property float ny\r\n";
    myfile << "property float nz\r\n";
    myfile << "element face " << inputFaces.size() << "\r\n";
    myfile << "property list uchar uint vertex_indices\r\n";
    myfile << "end_header\r\n";

    for(unsigned long  i=0; i < inputCloud.size(); i++ ){
        myfile << inputCloud[i](0) <<" "<< inputCloud[i](1) <<" " << inputCloud[i](2) << " " <<
                  inputNormals[i](0) <<" "<< inputNormals[i](1) <<" " << inputNormals[i](2) << "\r\n";
    }
    for(unsigned long  i=0; i < inputFaces.size(); i++ ){
        myfile << "4" <<" "<< inputFaces[i](0) <<" "<< inputFaces[i](1) <<" " << inputFaces[i](2)<<" "<< inputFaces[i](3) <<"\r\n";
    }
    myfile.close();
    return 0;

}
FirstLayer::~FirstLayer(){

}

