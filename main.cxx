#include <algorithm>
#include <iostream>

#include "itkImage.h"
#include "itkOrientImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

using ImageTypeFloat3D = itk::Image< float, 3 >;

/**
\brief Function that performs orientation fix and returns the original orientation of the input image
*/
template< class TImageType = ImageTypeFloat3D >
std::pair< std::string, typename TImageType::Pointer > GetImageOrientation(const typename TImageType::Pointer inputImage, const std::string &desiredOrientation = "RAI")
{
  if (TImageType::ImageDimension != 3)
  {
    std::cerr << "This function is only defined for 3D images.\n";
    exit(EXIT_FAILURE);
  }
  // Map between axis string labels and SpatialOrientation
  std::map< std::string, itk::SpatialOrientation::ValidCoordinateOrientationFlags > orientationMap;
  orientationMap["RIP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP;
  orientationMap["LIP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP;
  orientationMap["RSP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP;
  orientationMap["LSP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP;
  orientationMap["RIA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA;
  orientationMap["LIA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA;
  orientationMap["RSA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA;
  orientationMap["LSA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA;
  orientationMap["IRP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP;
  orientationMap["ILP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP;
  orientationMap["SRP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP;
  orientationMap["SLP"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP;
  orientationMap["IRA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA;
  orientationMap["ILA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA;
  orientationMap["SRA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA;
  orientationMap["SLA"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA;
  orientationMap["RPI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI;
  orientationMap["LPI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI;
  orientationMap["RAI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;
  orientationMap["LAI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI;
  orientationMap["RPS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS;
  orientationMap["LPS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS;
  orientationMap["RAS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
  orientationMap["LAS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS;
  orientationMap["PRI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI;
  orientationMap["PLI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI;
  orientationMap["ARI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI;
  orientationMap["ALI"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI;
  orientationMap["PRS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS;
  orientationMap["PLS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS;
  orientationMap["ARS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS;
  orientationMap["ALS"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS;
  orientationMap["IPR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR;
  orientationMap["SPR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR;
  orientationMap["IAR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR;
  orientationMap["SAR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR;
  orientationMap["IPL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL;
  orientationMap["SPL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL;
  orientationMap["IAL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL;
  orientationMap["SAL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL;
  orientationMap["PIR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR;
  orientationMap["PSR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR;
  orientationMap["AIR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR;
  orientationMap["ASR"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR;
  orientationMap["PIL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL;
  orientationMap["PSL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL;
  orientationMap["AIL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL;
  orientationMap["ASL"] = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL;

  auto orientFilter = itk::OrientImageFilter< TImageType, TImageType >::New();
  orientFilter->SetInput(inputImage);
  orientFilter->UseImageDirectionOn();
  orientFilter->SetDirectionTolerance(0);
  orientFilter->SetCoordinateTolerance(0);

  auto desiredOrientation_wrap = desiredOrientation;
  std::transform(desiredOrientation_wrap.begin(), desiredOrientation_wrap.end(), desiredOrientation_wrap.begin(), ::toupper);

  // set the desired orientation and update
  orientFilter->SetDesiredCoordinateOrientation(orientationMap[desiredOrientation_wrap]);
  orientFilter->Update();
  auto outputImage = orientFilter->GetOutput();

  std::string originalOrientation;

  auto originalOrientation_code = orientFilter->GetGivenCoordinateOrientation();
  for (auto it = orientationMap.begin(); it != orientationMap.end(); ++it)
  {
    if (it->second == originalOrientation_code)
    {
      originalOrientation = it->first;
      break;
    }
  }
  if (originalOrientation.empty())
  {
    originalOrientation = "Unknown";
  }

  return std::make_pair(originalOrientation, outputImage);
}

int main(int argc, char** argv)
{
  if (argc < 2 || argc > 3)
  {
    std::cout << "Usage:\n\n" <<
      "ITKOrientFilterExample ${inputImageFile_nifti} ${outputImageFile_nifti}\n";
    return EXIT_FAILURE;
  }

  auto inputImageFile = std::string(argv[1]);
  auto outputImageFile = std::string(argv[2]);

  auto reader_1 = itk::ImageFileReader< ImageTypeFloat3D >::New();
  reader_1->SetFileName(inputImageFile);
  try
  {
    reader_1->Update();
  }
  catch (const std::exception&e)
  {
    std::cerr << "Input Image wasn't read: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  
  auto inputImage = reader_1->GetOutput();
  std::pair< std::string, ImageTypeFloat3D::Pointer > output;
  try
  {
    output = GetImageOrientation(inputImage);
  }
  catch (const std::exception&e)
  {
    std::cerr << "Couldn't orient the image properly: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  std::cout << "Original Orientation: " << output.first << "\n";
  
  auto writer_1 = itk::ImageFileWriter< ImageTypeFloat3D >::New();
  writer_1->SetFileName(outputImageFile);
  writer_1->SetInput(output.second);
  try
  {
    writer_1->Update();
  }
  catch (const std::exception&e)
  {
    std::cerr << "Couldn't write the image: " << e.what() << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "Started verification.\n";
  auto reader_2 = itk::ImageFileReader< ImageTypeFloat3D >::New();
  reader_2->SetFileName(outputImageFile);
  reader_2->Update(); // we know that the image at this stage 'should' be valid

  auto inputImage_rai_verify = reader_2->GetOutput();

  auto output_verify_fromRaw = GetImageOrientation(output.second, output.first);
  auto output_verify_fromFile = GetImageOrientation(inputImage_rai_verify, output.first);

  auto inputImage_rai_verify_original_fromRaw = output_verify_fromRaw.second;
  auto inputImage_rai_verify_original_fromFile = output_verify_fromFile.second;

  std::cout << "Original Image Properties:\n";
  std::cout << "\tOrigin: " << inputImage->GetOrigin() << "\n" 
    //<< "\tDirection Cosines: " << inputImage->GetDirection() << "\n"
    ;

  std::cout << "Oriented Image [RAI] Properties:\n";
  std::cout << "\tOrigin: " << inputImage_rai_verify->GetOrigin() << "\n"
    //<< "\tDirection Cosines: " << inputImage_rai_verify->GetDirection() << "\n"
    ;

  std::cout << "Re-Oriented Image [LPI] - from Raw RAI Properties:\n";
  std::cout << "\tOrigin: " << inputImage_rai_verify_original_fromRaw->GetOrigin() << "\n"
    //<< "\tDirection Cosines: " << inputImage_rai_verify_original->GetDirection() << "\n"
    ;

  std::cout << "Re-Oriented Image [LPI] - from Written RAI Properties:\n";
  std::cout << "\tOrigin: " << inputImage_rai_verify_original_fromFile->GetOrigin() << "\n"
    //<< "\tDirection Cosines: " << inputImage_rai_verify_original->GetDirection() << "\n"
    ;

  std::cout << "Finished successfully.\n";
  return EXIT_SUCCESS;
}

/// code from mihail.isakov
//template<typename T> bool orient_filter(
//  const typename T::Pointer & image,
//  typename T::Pointer & out_image,
//  int &originalOrientation,
//  int f)
//{
//  if (image.IsNull()) return false;
//  typedef itk::OrientImageFilter<T, T> OrientImageFilterType;
//  typename OrientImageFilterType::Pointer filter = OrientImageFilterType::New();
//  try
//  {
//    filter->SetInput(image);
//    filter->UseImageDirectionOn();
//    filter->SetDesiredCoordinateOrientation(
//      static_cast<itk::SpatialOrientation::ValidCoordinateOrientationFlags>(f));
//    filter->Update();
//    out_image = filter->GetOutput();
//  }
//  catch (itk::ExceptionObject & ex)
//  {
//    std::cout << ex.GetDescription() << std::endl;
//    return false;
//  }
//  if (out_image.IsNotNull()) out_image->DisconnectPipeline();
//  else return false;
//  
//  originalOrientation = filter->GetGivenCoordinateOrientation();
//  return true;
//}
//
//int main(int argc, char ** argv)
//{
//  typedef itk::Image<float, 3> ImageType;
//  typedef itk::ImageFileReader<ImageType> ReaderType;
//  typedef itk::ImageFileWriter<ImageType> WriterType;
//  ImageType::Pointer inputLPI;
//  ImageType::Pointer outputRAI;
//  ImageType::Pointer outputLPI;
//  ImageType::Pointer reopenedRAI;
//  ImageType::Pointer reopenedLPI;
//  ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName(argv[1]);
//  try
//  {
//    reader->Update();
//    inputLPI = reader->GetOutput();
//  }
//  catch (itk::ExceptionObject & ex)
//  {
//    std::cout << ex.GetDescription() << std::endl;
//    return 1;
//  }
//  std::cout << "original LPI " << inputLPI->GetOrigin()
//    << std::endl;
//  int originalOrientation = 0, temp = 0;
//  orient_filter<ImageType>(
//    inputLPI,
//    outputRAI, originalOrientation,
//    (int)itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
//  std::cout << "filter output LPI->RAI " << outputRAI->GetOrigin()
//    << std::endl;
//  orient_filter<ImageType>(
//    outputRAI,
//    outputLPI, temp,
//    (int)originalOrientation);
//  std::cout << "filter output RAI->LPI " << outputLPI->GetOrigin()
//    << std::endl;
//  WriterType::Pointer writer1 = WriterType::New();
//  WriterType::Pointer writer2 = WriterType::New();
//  writer1->SetInput(outputRAI);
//  writer1->SetFileName("outputRAI.nii.gz");
//  writer2->SetInput(outputRAI);
//  writer2->SetFileName("outputLPI.nii.gz");
//  try
//  {
//    writer1->Update();
//    writer2->Update();
//  }
//  catch (itk::ExceptionObject & ex)
//  {
//    std::cout << ex.GetDescription() << std::endl;
//    return 1;
//  }
//  ReaderType::Pointer reader1 = ReaderType::New();
//  ReaderType::Pointer reader2 = ReaderType::New();
//  reader1->SetFileName("outputRAI.nii.gz");
//  reader2->SetFileName("outputLPI.nii.gz");
//  try
//  {
//    reader1->Update();
//    reader2->Update();
//    reopenedRAI = reader1->GetOutput();
//    reopenedLPI = reader2->GetOutput();
//  }
//  catch (itk::ExceptionObject & ex)
//  {
//    std::cout << ex.GetDescription() << std::endl;
//    return 1;
//  }
//  std::cout << "re-opened LPI->RAI " << reopenedRAI->GetOrigin()
//    << std::endl;
//  std::cout << "re-opened RAI->LPI n" << reopenedLPI->GetOrigin()
//    << std::endl;
//  return 0;
//}