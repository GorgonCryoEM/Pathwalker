/**
 * MRC data definition
 *
 * Author: Minxin Cheng
 * Date: 06/12/2012
 */


#pragma once

#include "general_data.h"

class MRC_Data :
  public General_Data
{
public:
  MRC_Data(View *_view3D);
  ~MRC_Data(void);
  void dataGeneration(char *filename);
  void buildGrid();

  void getScalarGP(int index, float *scalar);
  void getTensorGP(int index, float tensor[6]);
  void getGradientGP(int index, float gradient[3]);
  void getScalar(float position[3], float *scalar);
  void getTensor(float position[3], float tensor[6]);
  void getGradient(float position[3], float gradient[3]);
  void getScalarCubic(float position[3], float *scalar);
  void getTensorCubic(float position[3], float tensor[6]);
  void getGradientCubic(float position[3], float gradient[3]);
  void getTensorCubicTable(int axis, int i, int j, int k, int index, float tensor[6]);
  void getGradientCubicTable(int axis, int i, int j, int k, int index, float gradient[3]);
  void getGradientBicubicTable(int axis, int ii, int jj, int kk, float sGradients[subdNum][subdNum][3]);
private:
  void MRCGeneration(char *filename);
  void DTIGeneration(char *filename);

  friend class boost::serialization::access;
  template <class Archive> void save(Archive &ar, const unsigned int version) const
  {
    ar & boost::serialization::base_object<General_Data>(*this);
    if (resultDebug)
    {
      ar & sizex;
      ar & sizey;
      ar & sizez;
      for ( int i = 0 ; i < sizex * sizey * sizez ; i ++ )
      {
        ar & scalars[i];
      }
      for ( int i = 0 ; i < 6 * sizex * sizey * sizez ; i ++ )
      {
        ar & tensors[i];
      }
      for ( int i = 0 ; i < 3 * sizex * sizey * sizez ; i ++ )
      {
        ar & gradients[i];
      }
    }
  }

  template <class Archive> void load(Archive &ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<General_Data>(*this);
    if (resultDebug)
    {
      ar & sizex;
      ar & sizey;
      ar & sizez;
      scalars = new float [ sizex * sizey * sizez ];
      tensors = new float [ 6 * sizex * sizey * sizez ];
      gradients = new float [ 3 * sizex * sizey * sizez ];
      for ( int i = 0 ; i < sizex * sizey * sizez ; i ++ )
      {
        ar & scalars[i];
      }
      for ( int i = 0 ; i < 6 * sizex * sizey * sizez ; i ++ )
      {
        ar & tensors[i];
      }
      for ( int i = 0 ; i < 3 * sizex * sizey * sizez ; i ++ )
      {
        ar & gradients[i];
      }
    }
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
  float *scalars;
  float *tensors;
  float *gradients;
  int iteration;
  float *edgeTable;
  float *faceTable;
};
