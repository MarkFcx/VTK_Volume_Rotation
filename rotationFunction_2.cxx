/*=========================================================================

  Program:   Visualization Toolkit
  Module:    rotationFunction.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
// ./rotationFunction_2 -DICOM /home/cxfeng/Desktop/My_Work/DicomData/12
=========================================================================*/
// VTK includes
#include "vtkBoxWidget.h"
#include "vtkCamera.h"
#include "vtkCommand.h"
#include "vtkColorTransferFunction.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageData.h"
#include "vtkImageResample.h"
#include "vtkMetaImageReader.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPlanes.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkXMLImageDataReader.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include <vtkMath.h>
#include <vtkImageMathematics.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>

#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAxesActor.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkSphereSource.h>

#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkImageMathematics.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkCharArray.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkVolumeMapper.h>
#include <vtkImageMathematics.h>
#include <vtkActor.h>
#include <vtkImageCast.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkShortArray.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <fstream>

#include <vtkBMPWriter.h>
#include <iostream>
#include <time.h>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
void PrintUsage()
{
  cout << "Usage: " << endl;
  cout << endl;
  cout << "  FixedPointVolumeRayCastMapperCT <options>" << endl;
  cout << endl;
  cout << "where options may include: " << endl;
  cout << endl;
  cout << "  -DICOM <directory>" << endl;
  cout << "  -VTI <filename>" << endl;
  cout << "  -MHA <filename>" << endl;
  cout << "  -DependentComponents" << endl;
  cout << "  -Clip" << endl;
  cout << "  -MIP <window> <level>" << endl;
  cout << "  -CompositeRamp <window> <level>" << endl;
  cout << "  -CompositeShadeRamp <window> <level>" << endl;
  cout << "  -CT_Skin" << endl;
  cout << "  -CT_Bone" << endl;
  cout << "  -CT_Muscle" << endl;
  cout << "  -FrameRate <rate>" << endl;
  cout << "  -DataReduction <factor>" << endl;
  cout << endl;
  cout << "You must use either the -DICOM option to specify the directory where" << endl;
  cout << "the data is located or the -VTI or -MHA option to specify the path of a .vti file." << endl;
  cout << endl;
  cout << "By default, the program assumes that the file has independent components," << endl;
  cout << "use -DependentComponents to specify that the file has dependent components." << endl;
  cout << endl;
  cout << "Use the -Clip option to display a cube widget for clipping the volume." << endl;
  cout << "Use the -FrameRate option with a desired frame rate (in frames per second)" << endl;
  cout << "which will control the interactive rendering rate." << endl;
  cout << "Use the -DataReduction option with a reduction factor (greater than zero and" << endl;
  cout << "less than one) to reduce the data before rendering." << endl;
  cout << "Use one of the remaining options to specify the blend function" << endl;
  cout << "and transfer functions. The -MIP option utilizes a maximum intensity" << endl;
  cout << "projection method, while the others utilize compositing. The" << endl;
  cout << "-CompositeRamp option is unshaded compositing, while the other" << endl;
  cout << "compositing options employ shading." << endl;
  cout << endl;
  cout << "Note: MIP, CompositeRamp, CompositeShadeRamp, CT_Skin, CT_Bone," << endl;
  cout << "and CT_Muscle are appropriate for DICOM data. MIP, CompositeRamp," << endl;
  cout << "and RGB_Composite are appropriate for RGB data." << endl;
  cout << endl;
  cout << "Example: FixedPointVolumeRayCastMapperCT -DICOM CTNeck -MIP 4096 1024" << endl;
  cout << endl;
}

int main(int argc, char *argv[])
{
	// Parse the parameters

	int count = 1;
	char *dirname = NULL;
	char *fileName=0;
	double opacityWindow = 4096;
	double opacityLevel = 2048;
	bool independentComponents=true;

	const int size = 512*512*600;
	const int size390 = 512*512*390;
//
	string dirP = "/home/cxfeng/Desktop/My_Work/BinaryData/Binary2/leftVol.bin";

	unsigned char* U_charArray_P = new unsigned char[size]; 
	char *binArray_P = new char[size];

	ifstream binFile2 (dirP.c_str(), ios::in | ios::binary);

	binFile2.read( binArray_P, size );
	binFile2.close();
	U_charArray_P = reinterpret_cast<unsigned char*>(binArray_P);

	vtkSmartPointer<vtkUnsignedCharArray> vtkUcharArray_P = 
	vtkSmartPointer<vtkUnsignedCharArray>::New();
	vtkUcharArray_P->SetArray(U_charArray_P, size, 1);

	vtkSmartPointer<vtkImageData> PelImagedata= 
	vtkSmartPointer<vtkImageData>::New();
	PelImagedata->GetPointData()->SetScalars(vtkUcharArray_P);
	PelImagedata->SetDimensions(512, 512, 600);
	PelImagedata->SetSpacing(0.96875, 0.96875, 1);
	PelImagedata->Modified();
//
    string dirL = "/home/cxfeng/Desktop/My_Work/BinaryData/102M/rightMaskSeg2.bin";
	unsigned char* U_charArray_L = new unsigned char[size390]; 
	char *binArray_L = new char[size390];

	ifstream binFile3 (dirL.c_str(), ios::in | ios::binary);

	binFile3.read( binArray_L, size );
	binFile3.close();
	U_charArray_L = reinterpret_cast<unsigned char*>(binArray_L);

	vtkSmartPointer<vtkUnsignedCharArray> vtkUcharArray_L = 
	vtkSmartPointer<vtkUnsignedCharArray>::New();
	vtkUcharArray_L->SetArray(U_charArray_L, size, 1);

	vtkSmartPointer<vtkImageData> LeftImagedata = 
	vtkSmartPointer<vtkImageData>::New();
	LeftImagedata->SetDimensions(512, 512, 390);
	LeftImagedata->SetSpacing(0.96875, 0.96875, 1);
	LeftImagedata->GetPointData()->SetScalars(vtkUcharArray_L);
	LeftImagedata->Modified();
//
    // string dirR = "/home/cxfeng/Desktop/My_Work/BinaryData/102M/rightMaskSeg2.bin";
	string dirR = "/home/cxfeng/Desktop/My_Work/BinaryData/102M/rightMaskSeg2.bin";
	unsigned char* U_charArray_R = new unsigned char[size]; 
    //short * U_charArray_R = new short[size];
	char *binArray_R = new char[size];

	ifstream binFile4 (dirR.c_str(), ios::in | ios::binary);

	binFile4.read( binArray_R, size );
	binFile4.close();
	U_charArray_R = reinterpret_cast<unsigned char*>(binArray_R);
	
	vtkSmartPointer<vtkUnsignedCharArray> vtkUcharArray_R = 
		vtkSmartPointer<vtkUnsignedCharArray>::New();
	vtkUcharArray_R->SetArray(U_charArray_R, size, 1);

	vtkSmartPointer<vtkImageData> RightImagedata = 
		vtkSmartPointer<vtkImageData>::New();
	RightImagedata->SetDimensions(512, 512, 390);
	RightImagedata->SetSpacing(.96875,.96875,1);
	RightImagedata->GetPointData()->SetScalars(vtkUcharArray_R);
	RightImagedata->Modified();

	
	// int s = 512*512*600;
	// char *buffer = new char[s]{0};
	// for(int i = 0; i < 512*100; i++){
	// 	buffer[i] = 1;
	// }
	// ofstream myFile("testwriteBinary2.bin",ios::out | ios::binary);
	// if(!myFile){
	// 	cout<<"false!!!"<<endl;
	// }
	// myFile.write(buffer, s);
	// myFile.close();
	//cout <<"RightImagedata: "<<RightImagedata->GetDataDimension()<<"\n";

//
	while ( count < argc )
	{
		if ( !strcmp( argv[count], "-DICOM" ) )
		{
		  size_t size = strlen(argv[count+1])+1;
		  dirname = new char[size];
		  snprintf( dirname, size, "%s", argv[count+1] );
		  count += 2;
		}
		else
		{
		  cout << "Unrecognized option: " << argv[count] << endl;
		  cout << endl;
		  PrintUsage();
		  exit(EXIT_FAILURE);
		}
	}

	// Read the data
 	//vtkSmartPointer<vtkAlgorithm> reader=0;
  	vtkSmartPointer<vtkImageData> dicomImagedata1=0;
  	vtkSmartPointer<vtkImageData> dicomImagedata2=0;
    vtkSmartPointer<vtkImageData> dicomImagedata3=0;
	if(dirname)
	{
		vtkSmartPointer<vtkDICOMImageReader> dicomReader = vtkSmartPointer<vtkDICOMImageReader>::New();
	    dicomReader->SetDirectoryName(dirname);
	    //dicomReader->SetDataScalarTypeToUnsignedChar();
	    dicomReader->Update();
	    dicomImagedata1=dicomReader->GetOutput();
	    //
	    vtkSmartPointer<vtkDICOMImageReader> dicomReader2 = vtkSmartPointer<vtkDICOMImageReader>::New();
	    dicomReader2->SetDirectoryName(dirname);
	    //dicomReader->SetDataScalarTypeToUnsignedChar();
	    dicomReader2->Update();
	    dicomImagedata2=dicomReader2->GetOutput();

	    //
	    vtkSmartPointer<vtkDICOMImageReader> dicomReader3 = vtkSmartPointer<vtkDICOMImageReader>::New();
	    dicomReader3->SetDirectoryName(dirname);
	    //dicomReader->SetDataScalarTypeToUnsignedChar();
	    dicomReader3->Update();
	    dicomImagedata3=dicomReader3->GetOutput();

	    //reader=dicomReader;
	}

    int dim[3];
    dicomImagedata1->GetDimensions(dim);
    double spacing[3];
    dicomImagedata1->GetSpacing(spacing);

//	cout <<"RightImagedata:"<<endl;
//	RightImagedata->Print(cout);
//    cout <<"dicomImagedata:"<< endl;
//    dicomImagedata1->Print(cout);

/*
	vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
	cast->SetInputData(dicomImagedata);
	//cast->SetInputData(RightImagedata);
	cast->SetOutputScalarTypeToUnsignedChar();
	cast->Update();
	RightImagedata->Print(cout);
    dicomImagedata->Print(cout);
    cast->Print(cout);
*/
    /*


    dicomImagedata2->SetDimensions(512, 512, 600);
    dicomImagedata2->AllocateScalars(VTK_SHORT, 0);
    short *dcmPtr2 = (short*)dicomImagedata2->GetScalarPointer();
*/
    unsigned char *ptr1 = (unsigned char*)PelImagedata->GetScalarPointer();
    unsigned char *ptr2 = (unsigned char*)LeftImagedata->GetScalarPointer();
    unsigned char *ptr3 = (unsigned char*)RightImagedata->GetScalarPointer();
	short *dcmPtr1 = (short*)dicomImagedata1->GetScalarPointer();
	short *dcmPtr2 = (short*)dicomImagedata2->GetScalarPointer();
	short *dcmPtr3 = (short*)dicomImagedata3->GetScalarPointer();
//	
	/*
	clock_t t = clock();
	for(int i = 0; i < 512*512*600; i++)
	{

	    //*dcmPtr1 = (*dcmPtr1)*(*ptr1);
	    *dcmPtr2 = (*dcmPtr2)*(*ptr2);
	    *dcmPtr3 = (*dcmPtr3)*(*ptr3);

	    //dcmPtr1++;
	    dcmPtr2++;
	    dcmPtr3++;
	    
	    //ptr1++;	 
	    ptr2++;	 
	    ptr3++;	 

	}
	t = clock() - t;
	//dicomImagedata2->Modified();
 	cout << "Matrix* CPU time: " <<(double)t/CLOCKS_PER_SEC<<" s"<<"\n";
	*/
//

/*
	
	clock_t t = clock();
	for (int z = 0; z < 600; z++)
	{
		for (int y = 0; y < 512; y++)
		{
			for (int x = 0; x < 512; x++)
			{
			double* pixel = static_cast<double*>(dicomImagedata2->GetScalarPointer(x,y,z));
			*pixel = (*pixel)*(*(ptr3 + x + (512-y-1)*512 + z*512*512));
			//dcmPtr2++;
			}
		}
	}
	t = clock() - t;
	//dicomImagedata2->Modified();
 	cout << "Matrix* CPU time: " <<(double)t/CLOCKS_PER_SEC<<" s"<<"\n";

*/
    for (int z = 0; z < 600; z++)
    {
        for (int y = 0; y < 512; y++)
        {
            for (int x = 0; x < 512; x++)
            {
                if (!(*ptr1) == 0)
                {
                    //cout<<*ptr1<<endl;
                    *ptr1 = 2;
                    cout<<*ptr1<<endl;

                }
                ptr1++;
                //dcmPtr2++;
            }
        }
    }

//    for (int z = 0; z < 600; z++)
//    {
//        for (int y = 0; y < 512; y++)
//        {
//            for (int x = 0; x < 512; x++)
//            {
//                if (!(*ptr2) == 0)
//                {
//                    //cout<<*ptr2<<endl;
//                }
//                ptr2++;
//                //dcmPtr2++;
//            }
//        }
//    }
//
//	for (int z = 0; z < 600; z++)
//	{
//		for (int y = 0; y < 512; y++)
//		{
//			for (int x = 0; x < 512; x++)
//			{
//				if(!(*ptr1)==(*ptr2)){
//					cout<<"("<<x<<" "<<y<<" "<<z<<")  ";
//				}
//				ptr1++;
//				ptr2++;
//			//dcmPtr2++;
//			}
//			// cout<<endl;
//		}
//		// cout<<endl;
//	}
//
	/*
	clock_t t = clock();

	for(int i = 0; i < 600; i+=262144)
		for(int j = 511; j >= 0; --j)
			for(int z = 0; z < 512; ++z)
			{

			    //*dcmPtr1 = (*dcmPtr1)*(*ptr1);
			    *dcmPtr2 = (*dcmPtr2)*(*(ptr2+(512*j+z+i)));
			    *dcmPtr3 = (*dcmPtr3)*(*ptr3);

			    //dcmPtr1++;
			    dcmPtr2++;
			    dcmPtr3++;
			    
			    //ptr1++;	 
			    //ptr2++;	 
			    ptr3++;	 

			}
	t = clock() - t;
	//dicomImagedata2->Modified();
 	cout << "Matrix* CPU time: " <<(double)t/CLOCKS_PER_SEC<<" s"<<"\n";
	*/
//
 	/*
	for(int i = 0; i < 5; i++)
	{

	    //*dcmPtr1 = (*dcmPtr1)*(*ptr1);
	    *dcmPtr2 = (*dcmPtr2)*(*ptr2);
	    *dcmPtr3 = (*dcmPtr3)*(*ptr3);

	    //dcmPtr1++;
	    dcmPtr2++;
	    dcmPtr3++;
	    cout <<"dcmPtr: "<<dcmPtr2<<"\n";
	    
	    //ptr1++;	 
	    ptr2++;	 
	    ptr3++;
	    cout <<"ptr: "<<ptr2<<"\n";	 

	}
	*/

// Create our volume and mapper
//LeftView
	vtkSmartPointer<vtkVolume> volume11 = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkVolume> volume21 = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkVolume> volume31 = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkVolume> volume41 = vtkSmartPointer<vtkVolume>::New();

    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper11 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper21 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper31 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper41 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper11 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper21 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper31 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper41 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	
	mapper11->SetInputData(dicomImagedata1);
	mapper21->SetInputData(PelImagedata);
    mapper31->SetInputData(LeftImagedata);
    mapper41->SetInputData(RightImagedata);

//RightView
    vtkSmartPointer<vtkVolume> volume12 = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkVolume> volume22 = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkVolume> volume32 = vtkSmartPointer<vtkVolume>::New();
    vtkSmartPointer<vtkVolume> volume42 = vtkSmartPointer<vtkVolume>::New();

    vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper12 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper22 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper32 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> mapper42 = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper12 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper22 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper32 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> mapper42 = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	
	//
	mapper12->SetInputData(dicomImagedata1);
	mapper22->SetInputData(dicomImagedata1);
    mapper32->SetInputData(dicomImagedata2);
    mapper42->SetInputData(dicomImagedata3);

	mapper11->SetBlendModeToMaximumIntensity();
	mapper21->SetBlendModeToMaximumIntensity();
	mapper31->SetBlendModeToMaximumIntensity();
	mapper22->SetBlendModeToMaximumIntensity();
	mapper32->SetBlendModeToMaximumIntensity();
	mapper42->SetBlendModeToMaximumIntensity();

/*
    vtkSmartPointer<vtkImageMathematics> imagemath = vtkSmartPointer<vtkImageMathematics>::New();
    imagemath->SetOperationToMultiply();
    imagemath->SetInput1Data( RightImagedata );
    imagemath->SetInput2Data( dicomImagedata );
    imagemath->Update();
   
    int dims[3];
    imagemath->GetOutput()->GetDimensions(dims);//GetOutput()->
    cout <<"imagemath GetDimensions: "<< dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"\n";

    double spaceing[3];
    imagemath->GetOutput()->GetSpacing(spaceing);
    cout <<"imagemath GetSpaceing: "<< spaceing[0]<<" "<<spaceing[1]<<" "<<spaceing[2]<<"\n";
*/


//
    volume11->SetMapper( mapper11 );
    volume21->SetMapper( mapper21 );
    volume31->SetMapper( mapper31 );
    volume41->SetMapper( mapper41 );


    volume12->SetMapper( mapper12 );
    volume22->SetMapper( mapper22 );
    volume32->SetMapper( mapper32 );
    volume42->SetMapper( mapper42 );


//
	

//Rotation

	vtkSmartPointer<vtkTransform> trans1 = 
		vtkSmartPointer<vtkTransform>::New();
	trans1->PostMultiply();	
	clock_t t1 = clock();
    trans1->Translate(-155, -211.188, -83);
    trans1->RotateX(-8.19445);
    trans1->RotateY(100.0733765);
	trans1->RotateZ(3.75134);
	trans1->Translate(155, 211.188, 83);
    volume12->SetUserTransform(trans1);
    t1 = clock() - t1;

	vtkSmartPointer<vtkTransform> trans2 = 
		vtkSmartPointer<vtkTransform>::New();
	trans2->PostMultiply();
    clock_t t2 = clock();
    trans2->Translate(-155, -211.188, -83);
    trans2->RotateX(-8.19445);
    trans2->RotateY(0.0733765);
	trans2->RotateZ(3.75134);
	trans2->Translate(155, 211.188, 83);
    volume22->SetUserTransform(trans2);
    t2 = clock() - t2;
//    trans2->GetMatrix()->Print(cout);
    
	vtkSmartPointer<vtkTransform> trans3 = 
		vtkSmartPointer<vtkTransform>::New();
	trans3->PostMultiply();
	clock_t t3 = clock();
	trans3->Translate(200,200, 200);
    // trans3->Translate(-163.947, -255.421, -166.661);
 //    trans3->RotateX(10.5699);
 //    trans3->RotateY(0.484976);
	// trans3->RotateZ(18.9191);
	// trans3->Translate(163.947, 255.421, 166.661);
    volume32->SetUserTransform(trans3);
    t3 = clock() - t3;
//    trans3->GetMatrix()->Print(cout);
    
    vtkSmartPointer<vtkTransform> trans4 = 
		vtkSmartPointer<vtkTransform>::New();
	trans4->PostMultiply();
	clock_t t4 = clock();
    trans4->Translate(-340.088, -251.34, -157.613);
    trans4->RotateX(10.8356);
    trans4->RotateY(-0.788843);
	trans4->RotateZ(-7.88475);
	trans4->Translate(340.088, 251.34, 157.613);
	volume42->SetUserTransform(trans4);
	t4 = clock() - t4;
//	trans4->GetMatrix()->Print(cout);

	//cout << "DICOM     Rotation CPU time: " <<(double)t1/CLOCKS_PER_SEC<<"\n";
	cout << "Pelvis    Rotation CPU time: " <<(double)t2/CLOCKS_PER_SEC<<"\n";
	cout << "LeftMask  Rotation CPU time: " <<(double)t3/CLOCKS_PER_SEC<<"\n";
    cout << "RightMask Rotation CPU time: " <<(double)t4/CLOCKS_PER_SEC<<"\n";
    
//oriAxes and text
	#define LINE_LEN 100.
	vtkSmartPointer<vtkAxesActor> oriAxesActor =
		vtkSmartPointer<vtkAxesActor>::New();
	oriAxesActor->SetPosition(0, 0, 0);
	oriAxesActor->SetTotalLength(LINE_LEN, LINE_LEN, LINE_LEN);
	oriAxesActor->SetShaftType(0);
	oriAxesActor->SetAxisLabels(1);
	oriAxesActor->SetCylinderRadius(0.02);

	vtkSmartPointer<vtkTextActor> textActor1 = 
		vtkSmartPointer<vtkTextActor>::New();
	//textActor->SetPosition2(100, 40);
	textActor1->GetTextProperty()->SetFontSize(50);
	textActor1->GetTextProperty()->SetColor(0, 1, 1);

	vtkSmartPointer<vtkTextActor> textActor2 = 
		vtkSmartPointer<vtkTextActor>::New();
	//textActor->SetPosition2(100, 40);
	textActor2->GetTextProperty()->SetFontSize(50);
	textActor2->GetTextProperty()->SetColor(0, 1, 1);

	//textActor->SetInput("pelvis enter:\nx: 155.0 y: 211.188 z: 83.0\npelvis rotation angle:\nRotateX( -8.19445 )\nRotateY( 0.0733765 )\nRotateZ( 3.75134 )\n");
	textActor1->SetInput("Before rotation");
	textActor2->SetInput("After rotation");



//color and opacity and property
	// Create our transfer function
	vtkSmartPointer<vtkColorTransferFunction> colorFun = vtkSmartPointer<vtkColorTransferFunction>::New();
	vtkSmartPointer<vtkPiecewiseFunction> opacityFun = vtkSmartPointer<vtkPiecewiseFunction>::New();
	vtkSmartPointer<vtkVolumeProperty> property = vtkSmartPointer<vtkVolumeProperty>::New();
	
	property->SetIndependentComponents(independentComponents);
	property->SetColor( colorFun );
	property->SetScalarOpacity( opacityFun );
	property->SetInterpolationTypeToLinear();

	opacityFun->AddSegment(3,0,1000,1);
	colorFun->AddRGBSegment(0.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0 );

//
	vtkSmartPointer<vtkColorTransferFunction> colorFun2 = vtkSmartPointer<vtkColorTransferFunction>::New();
	vtkSmartPointer<vtkPiecewiseFunction> opacityFun2 = vtkSmartPointer<vtkPiecewiseFunction>::New();
	vtkSmartPointer<vtkVolumeProperty> property2 = vtkSmartPointer<vtkVolumeProperty>::New();
	
	property2->SetIndependentComponents(independentComponents);
	property2->SetColor( colorFun2 );
	property2->SetScalarOpacity( opacityFun2 );
	property2->SetInterpolationTypeToLinear();

	opacityFun2->AddSegment(0.999,0,2.0,1);
	colorFun2->AddRGBSegment(0.0, 1.0, .1, 1.0, 3.0, 1.0, 1.0, 1.0 );

//
	vtkSmartPointer<vtkColorTransferFunction> colorFun3 = vtkSmartPointer<vtkColorTransferFunction>::New();
	vtkSmartPointer<vtkPiecewiseFunction> opacityFun3 = vtkSmartPointer<vtkPiecewiseFunction>::New();
	vtkSmartPointer<vtkVolumeProperty> property3 = vtkSmartPointer<vtkVolumeProperty>::New();
	
	property3->SetIndependentComponents(independentComponents);
	property3->SetColor( colorFun3 );
	property3->SetScalarOpacity( opacityFun3 );
	property3->SetInterpolationTypeToLinear();

	opacityFun3->AddSegment(0.999,0,1.0000,1);
	colorFun3->AddRGBSegment(0.0, .0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0 );

//
	vtkSmartPointer<vtkColorTransferFunction> colorFun4 = vtkSmartPointer<vtkColorTransferFunction>::New();
	vtkSmartPointer<vtkPiecewiseFunction> opacityFun4 = vtkSmartPointer<vtkPiecewiseFunction>::New();
	vtkSmartPointer<vtkVolumeProperty> property4 = vtkSmartPointer<vtkVolumeProperty>::New();
	
	property4->SetIndependentComponents(independentComponents);
	property4->SetColor( colorFun4 );
	property4->SetScalarOpacity( opacityFun4 );
	property4->SetInterpolationTypeToLinear();

	opacityFun4->AddSegment(0.999,0,1.0000,1);
	colorFun4->AddRGBSegment(0.0, 1.0, 1.0, .0, 3.0, 1.0, 1.0, 1.0 );

//
    // connect up the volume to the property and the mapper
    volume11->SetProperty( property2 );
    volume21->SetProperty( property2 );
    volume31->SetProperty( property3 );
    volume41->SetProperty( property4 );

    volume12->SetProperty( property );
    volume22->SetProperty( property2 );
    volume32->SetProperty( property3 );
    volume42->SetProperty( property4 );

	
// Create the renderer, render window and interactor
	vtkSmartPointer<vtkRenderer> renderer1 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderer> renderer2 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> renderer3 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> renderer4 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer1);
	renWin->AddRenderer(renderer2);
    renWin->AddRenderer(renderer3);
    renWin->AddRenderer(renderer4);
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	//iren->SetDesiredUpdateRate(frameRate / (1+clip) );
	
	//iren->GetInteractorStyle()->SetDefaultRenderer(renderer1);
	//iren->GetInteractorStyle()->SetDefaultRenderer(renderer2);

	//double leftViewport[] = { 0.0, 0.0, 0.5, 1.0 };
	///double rightViewport[] = { 0.5, 0.0, 1.0, 1.0 };

	double leftUpViewport[] = { 0.0, 0.5, 0.5, 1.0 };
	double rightUpViewport[] = { 0.5, 0.5, 1.0, 1.0 };
	double leftDownViewport[] = { 0.0, 0.0, 0.5, 0.5 };
	double rightDownViewport[] = { 0.5, 0.0, 1.0, 0.5 };

	// Set the default window size
	renWin->SetSize(1000,1000);

    renderer1->AddActor(oriAxesActor);
    //renderer1->AddActor(sphereActor);
    //renderer1->AddVolume( volume11 );
    renderer1->AddVolume( volume21 );
    renderer1->AddVolume( volume31 );
    //renderer1->AddVolume( volume41 );
    //renderer1->AddActor2D(textActor1);
  
  
    renderer2->AddActor(oriAxesActor);
    //renderer2->AddActor(sphereActor);
    //renderer2->AddVolume( volume12 );
    //renderer2->AddVolume( volume22 );
    //renderer2->AddVolume( volume32 );
    //renderer2->AddVolume( volume42 );
    //renderer2->AddActor2D(textActor2);

    //renderer3->AddActor(volume33);
    

    renderer1->SetBackground(0,0,0);
    renderer2->SetBackground(0,0,0);
    renderer3->SetBackground(0,0,0);
    renderer4->SetBackground(0,0,0);
    renderer1->SetViewport(leftUpViewport);
    renderer2->SetViewport(rightUpViewport);
    renderer3->SetViewport(leftDownViewport);
    renderer4->SetViewport(rightDownViewport);
  
    //renderer1->ResetCamera();
    //renderer2->ResetCamera();
    //renderer3->ResetCamera();
    //renderer4->ResetCamera();

 
	renWin->Render();
	iren->Start();

	sleep(5);
    vtkSmartPointer<vtkWindowToImageFilter> filter = vtkSmartPointer<vtkWindowToImageFilter>::New();
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	filter->SetInput(renWin);
	writer->SetInputConnection(filter->GetOutputPort());
	writer->SetFileName("mask+mask");
	writer->SetFilePattern("png");
	writer->Write();
    
    // interact with data
    
    

    
    return 0;
}
