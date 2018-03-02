#ifndef PBC_H
#define PBC_H

#include "structure.hpp"

std::array<bool, 3> PBC_f{{false, false, false}};

//---------return nearest neighbor atom b-----------//
Vector3D Nearest_atom(const Vector3D& a, const Vector3D& b)
{
   Vector3D result(b);
   Vector3D center_b(b);
   int counter=0;
   
   const double max_errorx = Cell_length[0]*std::numeric_limits<double>::epsilon();
   const double max_errory = Cell_length[1]*std::numeric_limits<double>::epsilon();
   const double max_errorz = Cell_length[2]*std::numeric_limits<double>::epsilon();
   while((PBC_f[0] && (std::abs(a.x-result.x)>((Cell_length[0]+max_errorx)*0.5))) 
         || (PBC_f[1] && (std::abs(a.y-result.y)>((Cell_length[1]+max_errory)*0.5)))
         || (PBC_f[2] && (std::abs(a.z-result.z)>((Cell_length[2]+max_errorz)*0.5)))
         )
   {
      ++counter;
      for(int i=0; i<3; ++i)
      {
         for(int j=0; j<3; ++j)
         {
            for(int k=0; k<3; ++k)
            {
               Vector3D b_dush = center_b;
               if(PBC_f[0]){  if(i!=0){ (i==1)? b_dush.x +=Cell_length[0]: b_dush.x -=Cell_length[0]; } }
               if(PBC_f[1]){  if(j!=0){ (j==1)? b_dush.y +=Cell_length[1]: b_dush.y -=Cell_length[1]; } }
               if(PBC_f[2]){  if(k!=0){ (k==1)? b_dush.z +=Cell_length[2]: b_dush.z -=Cell_length[2]; } }
               if((a-result).norm() > (a-b_dush).norm())
               {
                  result = b_dush;
               }
            }
         }
      }
      center_b = result;
      if(counter>30)
      { 
         std::cout<<"Error! counter is over in PBC."<<std::endl;
         exit(1); 
      }
   }
   if(!PBC_f[0]){ result.x=b.x;}
   if(!PBC_f[1]){ result.y=b.y;}
   if(!PBC_f[2]){ result.z=b.z;}

   return result;
}

bool is_Box_in(const Vector3D& p)
{
   bool result=true;
   for(int i=0; i<3; ++i)
   {
      if(PBC_f[i])
      {
         if(p[i]<0 || p[i]>Cell_length[i]){  result=false;  }
      }
   }
   return result;
}

#endif
