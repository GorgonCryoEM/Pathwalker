/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by.edge
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#ifndef eman_processor_h__
#define eman_processor_h__ 1

#include "emobject.h"
#include "util.h"
#include "geometry.h"
#include "transform.h"
#include "emdata.h"
#include "gorgon/skeletonizer.h"

#include <cfloat>
#include <climits>
#include <cstring>

using std::vector;
using std::map;
using std::string;

namespace EMAN
{
	class EMData;

	/** Typical usage of Processors are as follows:
     *
     *   - How to get all the processor names
     *@code
     *      vector<string> all_processors = Factory<Processor>::get_list();
     @endcode
     *   - How to use a processor
     *@code
     *      EMData *img = ...;
     *      img->process_inplace("PROCESSORNAME", Dict("sigma", 12));
     @endcode
     *   - How to define a new XYZProcessor \n
     *      XYZProcessor should either extend the base class 'Processor' or a
     *      subclass of 'Processor'. At a minimum, it should define:
	 *      (Please replace 'XYZ' with your own class name).
	 *@code
     *          string get_name() const { return "processorname"; }
     *          static Processor *NEW() { return XYZProcessor(); }
	 @endcode
	 *      If XYZProcessor is a parent class, it should define:
	 *@code
	 *          static string get_group_desc();
	 @endcode
	 *      Otherwise, it should define:
	 *@code
	 *          string get_desc() const;
	 @endcode
     *      If XYZProcessor need parameters not defined by its parent
     *      class, it should define:
	 *@code
	 *          Dict get_params() const;
	 *          void set_params(const Dict & new_params);
     *          TypeDict get_param_types() const;
     @endcode
     */
	class Processor
	{
	  public:
		virtual ~ Processor()
		{
		}

		/** To process an image in-place.
		 * For those processors which can only be processed out-of-place, override this function
		 * to just print out some error message to remind user call the out-of-place version.
		 * @param image The image to be processed.
		 */
		virtual void process_inplace(EMData *image) = 0;

		/** To proccess an image out-of-place.
		 * For those processors which can only be processed out-of-place, override this function
		 * to give the right behavior.
		 * @param image The image will be copied, actual process happen on copy of image.
		 * @return the image processing result, may or may not be the same size of the input image
		 * */
		virtual EMData* process(const EMData * const image);

		/** To process multiple images using the same algorithm.
		 * @param images Multiple images to be processed.
		 */
		virtual void process_list_inplace(vector < EMData * > & images)
		{
			for (size_t i = 0; i < images.size(); i++) {
				process_inplace(images[i]);
			}
		}

		/** Get the processor's name. Each processor is identified by a unique name.
		 * @return The processor's name.
		 */
		virtual string get_name() const = 0;

		/** Get the processor parameters in a key/value dictionary.
		 * @return A key/value pair dictionary containing the parameters.
		 */
		virtual Dict get_params() const
		{
			return params;
		}

		/** Set the processor parameters using a key/value dictionary.
		 * @param new_params A dictionary containing the new parameters.
		 */
		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
		}

		/** Get processor parameter information in a dictionary. Each
		 * parameter has one record in the dictionary. Each record
		 * contains its name, data-type, and description.
		 *
		 * @return A dictionary containing the parameter info.
		 */
		virtual TypeDict get_param_types() const
		{
			return TypeDict();
		}

		/** Get the description of this group of processors. This
		 * function is defined in a parent class. It gives a
		 * introduction to a group of processors.
		 *
		 * @return The description of this group of processors.
		 */
		static string get_group_desc()
		{
			return "EMAN processors are in-place image processors. You may apply a processor to process a single image or process multiple images. Processor class is the base class for all processor. <br> \
The basic design of EMAN Processors: <br>\
    1) Each Processor class defines an image-processinging algorithm. <br>\
    2) All the Processor classes in EMAN are managed by a Factory pattern. So each Processor class must define: <br> a) a unique name to idenfity itself in the factory. <br>b) a static method to register itself in the factory.<br>\
    3) Each Processor class defines its own parameter set.<br>\
    4) Each Processor class defines functions to return its documentation including parameter information, and processor description. These functions enable EMAN to generate processor manuals dynamically.";
		}

		/** Get the descrition of this specific processor. This function
		 * must be overwritten by a subclass.
		 *
		 * @return The description of this processor.
		 */
		virtual string get_desc() const = 0;

		/** Fourier filter Processor type enum.
		 *  New Fourier filter processors are computed in a single function,
		 *  EMFourierFilterFunc, that uses a large switch statement to
		 *  apply the correct filter processor.  This enum specifies the
		 *  filter processor to be applied.
		 */
		enum fourier_filter_types {
			TOP_HAT_LOW_PASS,
			TOP_HAT_HIGH_PASS,
			TOP_HAT_BAND_PASS,
			TOP_HOMOMORPHIC,
			GAUSS_LOW_PASS,
			GAUSS_HIGH_PASS,
			GAUSS_BAND_PASS,
			GAUSS_INVERSE,
			GAUSS_HOMOMORPHIC,
			BUTTERWORTH_LOW_PASS,
			BUTTERWORTH_HIGH_PASS,
			BUTTERWORTH_HOMOMORPHIC,
			KAISER_I0,
			KAISER_SINH,
			KAISER_I0_INVERSE,
			KAISER_SINH_INVERSE,
			SHIFT,
			TANH_LOW_PASS,
			TANH_HIGH_PASS,
			TANH_HOMOMORPHIC,
			TANH_BAND_PASS,
			RADIAL_TABLE,
		        CTF_,
		};

		/** Compute a Fourier-filter processed image in place.
		 *
		 *  @par Purpose: Apply selected Fourier space processor to 1-,2-, or 3-D image.
		 *  @par Method:
		 *
		 *  @param     fimage  Input image object to be processed, either
		 *                     a real-space image or a Fourier-space image.
		 *                     Image may be 1-, 2-, or 3-dimensional.  The
		 *                     original input image is not touched by
		 *                     this routine.
		 *
		 *  @param[in] params  Processor parameters.  Different processors require
		 *                     different parameters, so we this routine accepts
		 *                     a dictionary of parameters and looks up the
		 *                     appropriate params for the chosen processor at
		 *                     run time.  All processors use the "dopad"
		 *                     parameter to determine whether the
		 *                     Fourier workspace array should be zero-
		 *                     padded to twice the original length
		 *                     (dopad == 1) or not zero-padded at all
		 *                     (dopad == 0).
		 *  @return No explicit return.  The image fimage is modified
		 *  in place.
		 */
		static void
		EMFourierFilterInPlace(EMData* fimage, Dict params) {
			bool doInPlace = true;
			EMFourierFilterFunc(fimage, params, doInPlace);
		}

		/** Compute a Fourier-processor processed image without altering the original image.
		 *
		 *  @par Purpose: Apply selected Fourier space processor to 1-,2-, or 3-D image.
		 *  @par Method:
		 *
		 *  @param     fimage  Input image object to be processeded, either
		 *                     a real-space image or a Fourier-space image.
		 *                     Image may be 1-, 2-, or 3-dimensional.
		 *
		 *  @param[in] params  Processor parameters.  Different processors require
		 *                     different parameters, so we this routine accepts
		 *                     a dictionary of parameters and looks up the
		 *                     appropriate params for the chosen processor processor at
		 *                     run time.  All processors use the "dopad"
		 *                     parameter to determine whether the
		 *                     Fourier workspace array should be zero-
		 *                     padded to twice the original length
		 *                     (dopad == 1) or not zero-padded at all
		 *                     (dopad == 0).
		 *  @return 1-, 2-, or 3-dimensional filter processed image.  If the
		 *          input image is a real-space image, then the returned
		 *          output image will also be a real-space image.
		 *          Similarly, if the input image is already a Fourier image,
		 *          then the output image will be a Fourier image.
		 */
		static EMData*
		EMFourierFilter(EMData* fimage, Dict params) {
			bool doInPlace = false;
			return EMFourierFilterFunc(fimage, params, doInPlace);
		}

	  private:
		/** Compute a Fourier-filter processed image.
		 *  This function is called by either of the convience functions
		 *  EMFourierFilter or EMFourierFilterInPlace.
		 *
		 *  @par Purpose: Apply selected Fourier space processor to 1-,2-, or 3-D image.
		 *  @par Method:
		 *
		 *  @param     fimage  Input image object to be processed, either
		 *                     a real-space image or a Fourier-space image.
		 *                     Image may be 1-, 2-, or 3-dimensional.  Image
		 *                     fimage will not be changed unless
		 *                     inplace == true.
		 *  @param[in] params  Processor parameters.  Different processor processors require
		 *                     different parameters, so we this routine accepts
		 *                     a dictionary of parameters and looks up the
		 *                     appropriate params for the chosen processor processor at
		 *                     run time.  All processors use the "dopad"
		 *                     parameter to determine whether the
		 *                     Fourier workspace array should be zero-
		 *                     padded to twice the original length
		 *                     (dopad == 1) or not zero-padded at all
		 *                     (dopad == 0).
		 *  @param[in] doInPlace Inplace flag.  If this flag is true then
		 *                     fimage will contain the processeded image
		 *                     when this function returns.
		 *
		 *  @return 1-, 2-, or 3-dimensional filter processed image.  If the
		 *          input image is a real-space image, then the returned
		 *          output image will also be a real-space image.
		 *          Similarly, if the input image is already a Fourier image,
		 *          then the output image will be a Fourier image.
		 *          In either case, if inplace == true then the output
		 *          image (pointer) will be the same as the input image
		 *          (pointer).
		 */
		static EMData*
		EMFourierFilterFunc(EMData* fimage, Dict params, bool doInPlace=true);

	  protected:
			mutable Dict params;
	};

	class ImageProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		static string get_group_desc()
		{
			return "An Image Processor defines a way to create a processor image. The processor image is used to multiply the input-image in the fourier space. ImageFilter class is the base class. Each specific ImageFilter class must define function create_processor_image(). ";
		}

	  protected:
		virtual EMData * create_processor_image() const = 0;
	};

	

	/** Segment a volume into ~n subvolumes using K-means classification
	 *
	 *@author Steve Ludtke
	 *@date 2008/11/03
	 *@param ctf[in] A Ctf object to use
	 */
	class KmeansSegmentProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		virtual EMData* process(const EMData * const image);
		void process_inplace( EMData * image);

		TypeDict get_param_types() const
		{
			TypeDict d ;
			d.put("nseg", EMObject::INT, "Number of segments to divide the image into. default=12" );
			d.put("thr",EMObject::FLOAT,"Isosurface threshold value. Pixels below this will not be segmented");
			d.put("ampweight",EMObject::INT,"If set, will weight centers by voxel amplitude. default = 1");
			d.put("maxsegsize",EMObject::FLOAT,"Maximum radial distance from segment center to member voxel. Default=10000");
			d.put("minsegsep",EMObject::FLOAT,"Minimum segment separation. Segments too close will trigger a reseed");
			d.put("maxiter",EMObject::FLOAT,"Maximum number of iterations to run before stopping. Default=100");
			d.put("maxvoxmove",EMObject::FLOAT,"Maximum number of voxels that can move before quitting. Default=25");
			d.put("verbose",EMObject::INT,"Be verbose while running");
			/**
			 *An option for pseudoatom generation in pathwalker. Instead of random seeding, seed on the gird initially.
			 *@author Muyuan Chen
			 *@date 2014/06/05
	         */
			d.put("pseudoatom",EMObject::BOOL,"Doing pseudoatom generation");
			d.put("sep",EMObject::FLOAT,"Separation distance, used only in pseudoatom generation. Default=3.78");
			return d;
		}

		static Processor *NEW()
		{
			return new KmeansSegmentProcessor();
		}

		string get_desc() const
		{
			return "Performs K-means segmentation on a volume. Note that this method uses random seeds, and thus will return different results each time it is run. Returned map contains number of segment for each voxel (or 0 for unsegmented voxels). Segmentation centers are stored in 'segmentcenters' attribute, consisting of a list of 3n floats in x,y,z triples.";
		}

		static const string NAME;

	};



	/**The base class for real space processor working on individual pixels. The processor won't consider the pixel's coordinates and neighbors.
	 */
	class RealPixelProcessor:public Processor
	{
	  public:
		RealPixelProcessor():value(0), maxval(1), mean(0), sigma(0)
		{
		}
		void process_inplace(EMData * image);

		virtual void set_params(const Dict & new_params)
		{
			params = new_params;
			if (params.size() == 1) {
				vector < EMObject > dict_values = params.values();
				value = dict_values[0];
			}
		}

		static string get_group_desc()
		{
			return "The base class for real space processor working on individual pixels. The processor won't consider the pixel's coordinates and neighbors.";
		}

	  protected:
		virtual void process_pixel(float *x) const = 0;
		virtual void calc_locals(EMData *)
		{
		}
		virtual void normalize(EMData *) const
		{
		}

		float value;
		float maxval;
		float mean;
		float sigma;
	};

	
	/**f(x) = x if x >= minval; f(x) = 0 if x < minval
	*@param minval
	 */
	class ToZeroProcessor:public RealPixelProcessor
	{
		public:
			string get_name() const
			{
				return NAME;
			}
			static Processor *NEW()
			{
				return new ToZeroProcessor();
			}
			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("minval", EMObject::FLOAT, "Everything below this value is set to zero");
				return d;
			}

			string get_desc() const
			{
				return "f(x) = x if x >= minval; f(x) = 0 if x < minval.";
			}

			static const string NAME;

		protected:
			inline void process_pixel(float *x) const
			{
				if (*x < value) {
					*x = 0;
				}
			}
	};

	/**f(x) = x if x <= maxval; f(x) = 0 if x > maxval
	 * @param maxval
	 */
	
	/**Base class for normalization processors. Each specific normalization processor needs to define how to calculate mean and how to calculate sigma.
	 */
	class NormalizeProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		static string get_group_desc()
		{
			return "Base class for normalization processors. Each specific normalization processor needs to define how to calculate mean and how to calculate sigma.";
		}

	  protected:
		virtual float calc_sigma(EMData * image) const;
		virtual float calc_mean(EMData * image) const = 0;
	};

	/**Normalize an image so its vector length is 1.0.
	 */
	class NormalizeUnitProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeUnitProcessor();
		}

		string get_desc() const
		{
			return "Normalize an image so its vector length is 1.0.";
		}

		static const string NAME;

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};

 	inline float NormalizeUnitProcessor::calc_mean(EMData *) const { return 0; }

	/**Normalize an image so its elements sum to 1.0 (fails if mean=0)
	 */
	class NormalizeUnitSumProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeUnitSumProcessor();
		}

		string get_desc() const
		{
			return "Normalize an image so its elements sum to 1.0 (fails if mean=0)";
		}

		static const string NAME;

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};

	inline float NormalizeUnitSumProcessor::calc_mean(EMData *) const { return 0; }


	/**do a standard normalization on an image.
	 */
	class NormalizeStdProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeStdProcessor();
		}

		string get_desc() const
		{
			return "do a standard normalization on an image.";
		}

		static const string NAME;

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**Uses a 1/0 mask defining a region to use for the zero-normalization.if no_sigma is 1, standard deviation not modified.
	 *@param mask the 1/0 mask defining a region to use for the zero-normalization
	 *@param no_sigma if this flag is zero, only average under the mask will be substracted. set this flag to 1, standard deviation not modified
	 */
	class NormalizeMaskProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		string get_desc() const
		{
			return "Uses a 1/0 mask defining a region to use for the zero-normalization.if no_sigma is 1, standard deviation not modified.";
		}

		static Processor *NEW()
		{
			return new NormalizeMaskProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("mask", EMObject::EMDATA, "the 1/0 mask defining a region to use for the zero-normalization");
			d.put("no_sigma", EMObject::INT, "if this flag is zero, only average under the mask will be substracted. set this flag to 1, standard deviation not modified");
			return d;
		}

		static const string NAME;

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};

	/**Normalize the image whilst also removing any ramps. Ramps are removed first, then mean and sigma becomes 0 and 1 respectively
	* This is essential Pawel Penczek's preferred method of particle normalization
	* @author David Woolford
	* @date Mid 2008
	*/
	class NormalizeRampNormVar: public Processor
	{
		public:
			string get_name() const
			{
				return NAME;
			}

			static Processor *NEW()
			{
				return new NormalizeRampNormVar();
			}

			string get_desc() const
			{
				return "First call filter.ramp on the image, then make the mean 0 and norm 1";
			}

			void process_inplace(EMData * image);

			static const string NAME;
	};

	/** Normalize the mass of the image assuming a density of 1.35 g/ml (0.81 Da/A^3).
	 * Only works for 3D images. Essentially a replica of Volume.C in EMAN1.
	 *@author David Woolford (a direct port of Steve Ludtke's code)
	 *@date 01/17/09
	 *@param apix Angstrom per pixel of the image. If not set will use the apix_x attribute of the image
	 *@param mass The approximate mass of protein/structure in kilodaltons
	 *@param thr The isosurface threshold which encapsulates the structure
	 */
	class NormalizeByMassProcessor: public Processor
	{
		public:
			string get_name() const
			{
				return NAME;
			}

			static Processor *NEW()
			{
				return new NormalizeByMassProcessor();
			}

			string get_desc() const
			{
				return "Normalize the mass of the image assuming a density of 1.35 g/ml (0.81 Da/A^3) (3D only)";
			}

			TypeDict get_param_types() const
			{
				TypeDict d;
				d.put("apix", EMObject::FLOAT,"Angstrom per pixel of the image. If not set will use the apix_x attribute of the image");
				d.put("mass", EMObject::FLOAT,"The approximate mass of protein/structure in kilodaltons");
				d.put("thr", EMObject::FLOAT,"The isosurface threshold which encapsulates the structure");
				d.put("verbose", EMObject::INT,"If set will give details about the normalization");
				return d;
			}

			void process_inplace(EMData * image);

			static const string NAME;
	};


	/**normalizes an image, mean value equals to edge mean.
	 */
	class NormalizeEdgeMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeEdgeMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, mean value equals to edge mean.";
		}

		static const string NAME;

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image, mean value equals to mean of 2 pixel circular border.
	 */
	class NormalizeCircleMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeCircleMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, mean value equals to mean of 2 pixel circular radius or of the circular border if no radius is set.";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("radius", EMObject::FLOAT,"Radius of 2 pixel circular border");
			return d;
		}

		static const string NAME;

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image, uses 2 pixels on left and right edge
	 */
	class NormalizeLREdgeMeanProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeLREdgeMeanProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image, uses 2 pixels on left and right edge";
		}

		static const string NAME;

	  protected:
		float calc_mean(EMData * image) const;
	};

	/**normalizes an image. mean -> (maxval-minval)/2; std dev = (maxval+minval)/2;
	 */
	class NormalizeMaxMinProcessor:public NormalizeProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeMaxMinProcessor();
		}

		string get_desc() const
		{
			return "normalizes an image. mean -> (maxval-minval)/2; std dev = (maxval+minval)/2;";
		}

		static const string NAME;

	  protected:
		float calc_sigma(EMData * image) const;
		float calc_mean(EMData * image) const;
	};

	/**normalizes each row in the image individually
	 */
	class NormalizeRowProcessor:public Processor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeRowProcessor();
		}

		string get_desc() const
		{
			return "normalizes each row in the image individually";
		}

		static const string NAME;

		void process_inplace(EMData * image);
	};

	/**Sorry for the pun. This processor will take a second image and try to filter/scale it to optimally subtract it
	from the original image. The idea here is that if you have an image with noise plus a linear-filter modified projection,
	that a good measure of the similarity of the image to the projection would be to try and remove the projection from
	the image as optimally as possible, then compute the standard deviation of what's left.

	Now you might say that if the total energy in the noisy image is normalized then this should be equivalent to just
	integrating the FSC, which is what we use to do the optimal subtraction in the first place. This would be true, but
	this "optimal subtraction" has other purposes as well, such as the e2extractsubparticles program.
	 * @param ref Reference image to subtract
	 * @param return_radial Will return the radial filter function applied to ref as filter_curve
	 */
	class SubtractOptProcessor:public Processor
	{
	  public:
		virtual void process_inplace(EMData *image);
		virtual EMData* process(const EMData * const image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new SubtractOptProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("ref", EMObject::EMDATA, "Reference image to subtract");
			d.put("actual", EMObject::EMDATA, "If specified, ref is used for normalization, but actual is subtracted.");
			d.put("low_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] low cut-off frequency.");
			d.put("high_cutoff_frequency", EMObject::FLOAT, "Absolute [0,0.5] high cut-off frequency.");
			d.put("ctfweight",EMObject::BOOL, "Filter the image by CTF before subtraction");
			d.put("return_fft",EMObject::BOOL, "Skips the final IFT, and returns the FFT of the subtracted image");
			d.put("return_subim", EMObject::BOOL, "Instead of returning the image after subtraction, returns the filtered image which would have been subtracted from the image.");
			d.put("return_radial", EMObject::BOOL, "Return the radial filter function as an attribute (filter_curve)");
			d.put("return_presigma", EMObject::BOOL, "Return the sigma of the pre-subtracted image in real-space with the specified filter applied as sigma_presub. This is an expensive option.");
			return d;
		}

		string get_desc() const
		{
			return "This will filter/scale 'ref' optimally and subtract it from image using ring dot products in Fourier space for normalization. Cutoff frequencies apply a bandpass tophat filter to the output.";
		}

		static const string NAME;
	};


	/**use least square method to normalize
	 * @param to reference image normalize to
	 * @param low_threshold only take into account the reference image's pixel value between high and low threshold (zero is ignored)
	 * @param high_threshold only take into account the reference image's pixel value between high and low threshold (zero is ignored)
	 */
	class NormalizeToLeastSquareProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new NormalizeToLeastSquareProcessor();
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("to", EMObject::EMDATA, "reference image normalize to");
			d.put("ignore_zero", EMObject::BOOL, "If set, ignores any pixels which are exactly zero in either image. Defaut = True.");
			d.put("ignore_lowsig", EMObject::FLOAT, "If >0, then any pixels closer to the mean than val*sigma in either image excluded");
			d.put("low_threshold", EMObject::FLOAT, "only take into account the reference image's pixel value between high and low threshold (zero is always ignored)");
			d.put("high_threshold", EMObject::FLOAT, "only take into account the reference image's pixel value between high and low threshold (zero is always ignored)");
			return d;
		}

		string get_desc() const
		{
			return "use least square method to normalize";
		}

		static const string NAME;
	};



	class CutToZeroProcessor:public RealPixelProcessor
	{
	  public:
		string get_name() const
		{
			return NAME;
		}
		static Processor *NEW()
		{
			return new CutToZeroProcessor();
		}
		TypeDict get_param_types() const
		{
			TypeDict d;
			d.put("minval", EMObject::FLOAT, "the value that will be set to zero - all values below will also be set to zero. Values above get minval subtracted from them" );
			return d;
		}

		string get_desc() const
		{
			return "f(x) = x-minval if x >= minval; f(x) = 0 if x < minval.";
		}

		static const string NAME;

	  protected:
		void process_pixel(float *x) const
		{
		        *x = *x - value;
			if (*x < 0) {
				*x = 0;
			}
		}
	};


#ifdef SPARX_USING_CUDA
	/* class MPI CUDA kmeans processor
	 * 2009-02-13 17:34:45 JB first version
	 * 2009-09-02 11:19:10 JB for MPI version
	 * python wrap for GPU cluster
	 */
	class MPICUDA_kmeans {
	public:
		MPICUDA_kmeans();
		~MPICUDA_kmeans();
		int setup(int extm, int extN, int extn, int extK, int extn_start);
		void append_flat_image(EMData* im, int pos);
		int init_mem(int numdev);
		void compute_im2();
		int random_ASG(long int rnd);
		vector<int> get_ASG();
		vector<int> get_asg();
		void compute_NC();
		vector<int> get_NC();
		void set_ASG(const vector <int>& ASG);
		void set_NC(const vector <int>& NC);
		int get_ct_im_mv();
		void set_T(float extT);
		float get_T();
		void compute_AVE();
		void set_AVE(EMData* im, int pos);
		vector<EMData*> get_AVE();
		int one_iter();
		//int one_iter_SSE();
		//int AVE_to_host();
		int one_iter_SA();
		vector<float> compute_ji();
		vector<float> compute_criterion(const vector <float>& Ji);
		int shutdown();
	private:
		// params
		int m;
		int N;
		int n;
		int K;
		int nb_part;
		int n_start;
		int size_im;
		int size_IM;
		int size_AVE;
		int size_dist;
		int BLOCK_SIZE;
		int NB;
		int ins_BLOCK;
		int ite;
		float T;
		// debug
		int ct_im_mv;
		// host memory
		float* h_IM;
		float* h_im;
		float* h_AVE;
		float* h_dist;
		float* h_AVE2;
		float* h_im2;
		unsigned short int* h_ASG;
		unsigned short int* h_asg;
		unsigned int* h_NC;
		int* params;
		float ttt;
		// device memory
		float* d_im;
		float* d_AVE;
		float* d_dist;
		//int init_dist(); // intial h_dist and d_dist for SSE
                float compute_tt();
	};

#endif //EMAN2_USING_CUDA

#if 0

	class XYZProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new XYZProcessor();
		}

		string get_desc() const
		{
			return "N/A";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};


#endif


#if 0

	class XYZProcessor:public Processor
	{
	  public:
		void process_inplace(EMData * image);

		string get_name() const
		{
			return NAME;
		}

		static Processor *NEW()
		{
			return new XYZProcessor();
		}

		string get_desc() const
		{
			return "N/A";
		}

		TypeDict get_param_types() const
		{
			TypeDict d;
			return d;
		}

		static const string NAME;
	};


#endif


	int multi_processors(EMData * image, vector < string > processornames);
	void dump_processors();
	map<string, vector<string> > dump_processors_list();
	map<string, vector<string> > group_processors();

	template <> Factory < Processor >::Factory();
}

#endif	//eman_filter_h__

/* vim: set ts=4 noet: */

