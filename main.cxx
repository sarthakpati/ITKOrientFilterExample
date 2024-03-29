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

  auto desiredOrientation_wrap = desiredOrientation;
  std::transform(desiredOrientation_wrap.begin(), desiredOrientation_wrap.end(), desiredOrientation_wrap.begin(), ::toupper);

  // set the desired orientation and update
  auto orientFilter = itk::OrientImageFilter< TImageType, TImageType >::New();
  orientFilter->SetInput(inputImage);
  orientFilter->UseImageDirectionOn();
  //orientFilter->SetDirectionTolerance(0);
  //orientFilter->SetCoordinateTolerance(0);
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

///// code from mihail.isakov
//template<typename T> 
//bool orient_filter(
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

  std::string ext;
  //ext = ".img"; // doesn't work
  //ext = ".nii.gz"; // doesn't work
  //ext = ".vtk"; // doesn't work
  //ext = ".mha"; // works
  //ext = ".nrrd"; // works
  ext = ".mhd.gz"; // current test
  outputImageFile = "rai"; 
  auto writer_1 = itk::ImageFileWriter< ImageTypeFloat3D >::New();
  writer_1->SetFileName(outputImageFile + ext);
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
  reader_2->SetFileName(outputImageFile + ext);
  reader_2->Update(); // we know that the image at this stage 'should' be valid

  auto inputImage_rai_verify = reader_2->GetOutput();

  auto output_verify_fromRaw = GetImageOrientation(output.second, output.first);
  auto output_verify_fromFile = GetImageOrientation(inputImage_rai_verify, output.first);

  auto inputImage_rai_verify_original_fromRaw = output_verify_fromRaw.second;
  auto inputImage_rai_verify_original_fromFile = output_verify_fromFile.second;

  outputImageFile = "rai-lpi";
  auto writer_2 = itk::ImageFileWriter< ImageTypeFloat3D >::New();
  writer_2->SetFileName(outputImageFile + ext);
  writer_2->SetInput(inputImage_rai_verify_original_fromFile);
  writer_2->Update();

  std::cout << "Original Image Properties:\n";
  std::cout << "\tOrigin: " << inputImage->GetOrigin() << "\n" 
    //<< "\tDirection Cosines: " << inputImage->GetDirection() << "\n"
    ;

  std::cout << "Oriented Image [RAI] Raw Properties:\n";
  std::cout << "\tOrigin: " << output.second->GetOrigin() << "\n"
    //<< "\tDirection Cosines: \n" << output.second->GetDirection() << "\n"
    ;

  std::cout << "Oriented Image [RAI] FromFile Properties:\n";
  std::cout << "\tOrigin: " << inputImage_rai_verify->GetOrigin() << "\n"
    //<< "\tDirection Cosines:\n" << inputImage_rai_verify->GetDirection() << "\n"
    ;

  std::cout << "Re-Oriented Image [LPI] - from Raw RAI Properties:\n";
  std::cout << "\tOrigin: " << inputImage_rai_verify_original_fromRaw->GetOrigin() << "\n"
    //<< "\tDirection Cosines: " << inputImage_rai_verify_original->GetDirection() << "\n"
    ;

  std::cout << "Re-Oriented Image [LPI] - from FromFile RAI Properties:\n";
  std::cout << "\tOrigin: " << inputImage_rai_verify_original_fromFile->GetOrigin() << "\n"
    //<< "\tDirection Cosines: " << inputImage_rai_verify_original->GetDirection() << "\n"
    ;

  ///// testing mihail's code
  //ImageTypeFloat3D::Pointer orientedImage_rai_m, orientedImage_rai_lpi_m, orientedImage_rai_lpi_fromFile_m;
  //int originalOrintation_m, shouldGetRAI;
  //auto temp1 = orient_filter< ImageTypeFloat3D >(inputImage, orientedImage_rai_m, originalOrintation_m, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
  //auto temp2 = orient_filter< ImageTypeFloat3D >(orientedImage_rai_m, orientedImage_rai_lpi_m, shouldGetRAI, originalOrintation_m);

  //auto writer_2 = itk::ImageFileWriter< ImageTypeFloat3D >::New();
  //writer_2->SetFileName("test_m.nii.gz");
  //writer_2->SetInput(orientedImage_rai_m);
  //try
  //{
  //  writer_2->Update();
  //}
  //catch (const std::exception&e)
  //{
  //  std::cerr << "Couldn't write the image: " << e.what() << "\n";
  //  return EXIT_FAILURE;
  //}

  //auto reader_3 = itk::ImageFileReader< ImageTypeFloat3D >::New();
  //reader_3->SetFileName("test_m.nii.gz");
  //reader_3->Update(); // we know that the image at this stage 'should' be valid

  //auto temp3 = orient_filter< ImageTypeFloat3D >(reader_3->GetOutput(), orientedImage_rai_lpi_fromFile_m, shouldGetRAI, originalOrintation_m);

  //std::cout << "[M] Original Image Properties:\n";
  //std::cout << "\tOrigin: " << inputImage->GetOrigin() << "\n"
  //  //<< "\tDirection Cosines: " << inputImage->GetDirection() << "\n"
  //  ;

  //std::cout << "[M] Oriented Image [RAI] Raw Properties:\n";
  //std::cout << "\tOrigin: " << orientedImage_rai_m->GetOrigin() << "\n"
  //  //<< "\tDirection Cosines: \n" << orientedImage_rai_m->GetDirection() << "\n"
  //  ;

  //std::cout << "[M] Oriented Image [RAI] FromFile Properties:\n";
  //std::cout << "\tOrigin: " << inputImage_rai_verify->GetOrigin() << "\n"
  //  //<< "\tDirection Cosines:\n" << inputImage_rai_verify->GetDirection() << "\n"
  //  ;

  //std::cout << "[M] Re-Oriented Image [LPI] - from Raw RAI Properties:\n";
  //std::cout << "\tOrigin: " << reader_3->GetOutput()->GetOrigin() << "\n"
  //  //<< "\tDirection Cosines: " << inputImage_rai_verify_original->GetDirection() << "\n"
  //  ;

  //std::cout << "[M] Re-Oriented Image [LPI] - from FromFile RAI Properties:\n";
  //std::cout << "\tOrigin: " << orientedImage_rai_lpi_fromFile_m->GetOrigin() << "\n"
  //  //<< "\tDirection Cosines: " << inputImage_rai_verify_original->GetDirection() << "\n"
  //  ;

  ///// testing mihail's code

  std::cout << "Finished successfully.\n";
  return EXIT_SUCCESS;
}
