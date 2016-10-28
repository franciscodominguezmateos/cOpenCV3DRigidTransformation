/*
 * rigid_transformation.h
 *
 *  Created on: Jun 19, 2016
 *      Author: Francisco Dominguez
 *      A little library to work out exact rigid body transformation from paired 3d points
 */

#ifndef RIGID_TRANSFORMATION_H_
#define RIGID_TRANSFORMATION_H_

#include "opencv2/core/core.hpp"

using namespace std;
using namespace cv;

Mat centroidR(Mat &m){
	Mat t;
	m.row(0).copyTo(t);
	int n=m.rows;
	for(int i=1;i<n;i++){
		t+=m.row(i);
	}
	t/=n;
	return t;
}
Mat centroidC(Mat &m){
	Mat t=m.col(0);
	int n=m.cols;
	for(int i=1;i<n;i++){
		t+=m.col(i);
	}
	t/=n;
	return t;
}
Mat centerR(Mat &M,Mat &centroidM){
	Mat Ret;
	M.copyTo(Ret);
	int n=M.rows;
	for(int i=0;i<n;i++){
		Ret.row(i)-=centroidM;
	}
	return Ret;
}
Mat centerC(Mat &M){
	Mat Ret;
	M.copyTo(Ret);
	Mat centroidM(centroidC(M));
	int n=M.cols;
	for(int i=0;i<n;i++){
		Ret.col(i)-=centroidM;
	}
	return Ret;
}
void rigidTransformation(Mat &A,Mat &B,Mat &R,Mat &t){
	// from http://nghiaho.com/?page_id=671
	Mat centroidA,centroidB;
	centroidA=centroidR(A);
	centroidB=centroidR(B);
	Mat AA=centerR(A,centroidA);
	Mat BB=centerR(B,centroidB);
	Mat H=AA.t()*BB;
	Mat w, u, vt;
	SVD::compute(H, w, u, vt);
	R=vt.t()*u.t();
	if(cv::determinant(R)<0.0){
		vt.row(2)*=-1;
		R=vt.t()*u.t();
	}
	t=-R*centroidA.t()+centroidB.t();
}

#endif /* RIGID_TRANSFORMATION_H_ */
