import ants
import shutil
import argparse

class AntsRegistration:
    def __init__(self, fixed_image_path, moving_image_path=None, output_prefix=None, transform_files=None):
        self.fixed_image_path = fixed_image_path
        self.moving_image_path = moving_image_path
        self.output_prefix = output_prefix
        self.transform_files = transform_files
        self.transform = None

    def register(self, clean = False):
        if self.moving_image_path is None:
            raise ValueError("Moving image path must be provided for registration.")
        
        fixed_image = ants.image_read(self.fixed_image_path)
        moving_image = ants.image_read(self.moving_image_path)
        self.transform = ants.registration(fixed=fixed_image, moving=moving_image, type_of_transform='SyN')

        if not clean:
            ants.image_write(self.transform['warpedmovout'], f"{self.output_prefix}_warped.nii.gz")
            shutil.copy(self.transform['fwdtransforms'][0], f"{self.output_prefix}_0GenericAffine.mat")
            shutil.copy(self.transform['fwdtransforms'][1], f"{self.output_prefix}_1Warp.nii.gz")

    def apply_transform(self, image_path, output_path, is_label=False):
        if self.transform:
            transformlist = self.transform['fwdtransforms']
        elif self.transform_files:
            transformlist = self.transform_files
        else:
            raise ValueError("You need to run the registration first or provide transform files.")
        
        fixed_image = ants.image_read(self.fixed_image_path)
        image = ants.image_read(image_path)
        
        interpolation = 'nearestNeighbor' if is_label else 'linear'
        transformed_image = ants.apply_transforms(fixed=fixed_image, moving=image, transformlist=transformlist, interpolator=interpolation)
        #if is_label:
        #    transformed_image = transformed_image.astype('uint8')
        ants.image_write(transformed_image, output_path)

def main():
    parser = argparse.ArgumentParser(description="Perform non-linear registration using ANTsPy.")
    parser.add_argument('--fixed_image', type=str, help='Path to the fixed image')
    parser.add_argument('--moving_image', type=str, help='Path to the moving image')
    parser.add_argument('--output_prefix', type=str, help='Prefix for the output files')
    parser.add_argument('--apply_transform', type=str, help='Path to the image to apply the transform to')
    parser.add_argument('--output_transformed', type=str, help='Path to save the transformed image')
    parser.add_argument('--transform_files', nargs=2, help='Paths to the transform files (_0GenericAffine.mat and _1Warp.nii.gz)')
    parser.add_argument('--is_label', action='store_true', help='Indicate if the image to apply the transform to is a label mask')
    args = parser.parse_args()

    registration = AntsRegistration(
        fixed_image_path=args.fixed_image,
        moving_image_path=args.moving_image,
        output_prefix=args.output_prefix,
        transform_files=args.transform_files
    )

    if args.moving_image and args.output_prefix:
        registration.register()
    elif args.moving_image:
        registration.register(clean=True)

    if args.apply_transform and args.output_transformed:
        registration.apply_transform(args.apply_transform, args.output_transformed, args.is_label)

if __name__ == "__main__":
    main()