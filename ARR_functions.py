import SimpleITK as sitk
import os

def resampling(path_fixed_image, path_moving_image, output_folder, name_moving_image):
    fixed_image = sitk.ReadImage(path_fixed_image, sitk.sitkFloat32)
    moving_image = sitk.ReadImage(path_moving_image, sitk.sitkFloat32)
    output_dir = output_folder
    dimension = 3
    identity = sitk.Transform(dimension, sitk.sitkIdentity)
    moving_resampled = sitk.Resample(moving_image, fixed_image, identity, sitk.sitkNearestNeighbor, 0.0, moving_image.GetPixelID())
    sitk.WriteImage(moving_resampled, os.path.join(output_dir,  name_moving_image + '_resampled.nrrd'))
