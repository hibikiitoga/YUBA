#ifndef OUTPUT
#define OUTPUT
#include "Quaternion.hpp"

void Beads_base(   const unsigned int& step,const std::vector<Coordinate>& beads,const std::vector<double>& scalars);
void Beads_ex_base(const unsigned int& step,const std::vector<std::pair<Vector3D,double> >& beads_ex,const std::vector<double>& scalars);

void Beads(   const unsigned int& step,const std::vector<std::string>& list);
void Beads_ex(const unsigned int& step,const std::vector<std::string>& list);


void Membrane_base
(
   const unsigned int& step,
   const unsigned int& shf_mode,
   const std::vector<Vector3D>& vertices,
   const std::vector<Triangle>& triangles,
   const std::vector<double>& scalars
);

std::string                 input_file;
std::string                output_file;
double                       Bead_SIZE;
std::tuple<std::vector<Vector3D>,std::vector<Triangle>,std::vector<Quadrangle> > Sphere_template_ex;

void Beads_base
(
   const unsigned int& step,
   const std::vector<Coordinate>& beads,
   const std::vector<double>& scalars
)
{
   const std::string str_step = boost::lexical_cast<std::string> (step);
   const std::vector<Vector3D>& points = std::get<0>(Sphere_template_ex);
   const std::vector<Triangle>& triangles = std::get<1>(Sphere_template_ex);
   const std::vector<Quadrangle>& quadrangles = std::get<2>(Sphere_template_ex);
   std::string beads_file = output_file+"_beads_"+str_step+".vtk";
   std::ofstream ofs_b(beads_file,std::ios::trunc);
   ofs_b<<"# vtk DataFile Version 2.0"<<std::endl;
   ofs_b<<"Beads "<<input_file<<std::endl;
   ofs_b<<"ASCII"<<std::endl;
   ofs_b<<"DATASET POLYDATA"<<std::endl;
   ofs_b<<"POINTS"<<" "<<points.size()*beads.size()<<" "<<"float"<<std::endl;
   boost::format point("%f %f %f");
   for(size_t i=0,size=beads.size();i<size;++i)
   {
      const Vector3D& center = beads[i].v;
      for(size_t s=0,ssze=points.size();s<ssze;++s)
      {
         const Vector3D p = (points[s]*Bead_SIZE*0.5+center);
         ofs_b<<point %(float)p.x %(float)p.y %(float)p.z<<std::endl;
      }
   }
   boost::format pole("3 %d %d %d");
   boost::format equator("4 %d %d %d %d");
   ofs_b<<"POLYGONS "<<(triangles.size()+quadrangles.size())*beads.size()<<" "<<(4*triangles.size()+5*quadrangles.size())*beads.size()<<std::endl;

   const int polygon_size = (points.size());
   for(size_t i=0,size=beads.size();i<size;++i)
   {
      for(size_t t=0,stze=triangles.size();t<stze;++t)
      {
         ofs_b<<pole %(triangles[t].a+i*polygon_size) %(triangles[t].b+i*polygon_size) %(triangles[t].c+i*polygon_size)<<std::endl;
      }
      for(size_t q=0,sqze=quadrangles.size();q<sqze;++q)
      {
         ofs_b<<equator %(quadrangles[q].a+i*polygon_size) %(quadrangles[q].b+i*polygon_size) %(quadrangles[q].c+i*polygon_size) %(quadrangles[q].d+i*polygon_size)<<std::endl;
      }
   }

   if(!scalars.empty())
   {
      ofs_b<<"POINT_DATA "<<points.size()*beads.size()<<std::endl;
      ofs_b<<"SCALARS "<<"beads "<<"float 1"<<std::endl;
      ofs_b<<"LOOKUP_TABLE "<<"default"<<std::endl;
      boost::format cl("%f ");
      for(size_t i=0, size=beads.size(); i<size; ++i)
      {
         for(size_t j=0, sjze=points.size(); j<sjze; ++j)
         {
            ofs_b<<cl %(float)scalars[i];
         }
         ofs_b<<std::endl;
      }
   }
   ofs_b.close();

   return ;
}

void Beads_ex_base
(
   const unsigned int& step,
   const std::vector<std::pair<Vector3D,double> >& beads_ex,
   const std::vector<double>& scalars
)
{
   const std::string str_step = boost::lexical_cast<std::string> (step);
   const std::vector<Vector3D>& points = std::get<0>(Sphere_template_ex);
   const std::vector<Triangle>& triangles = std::get<1>(Sphere_template_ex);
   const std::vector<Quadrangle>& quadrangles = std::get<2>(Sphere_template_ex);
   std::string beads_file = output_file+"_beads_"+str_step+".vtk";
   std::ofstream ofs_b(beads_file,std::ios::trunc);
   ofs_b<<"# vtk DataFile Version 2.0"<<std::endl;
   ofs_b<<"Beads "<<input_file<<std::endl;
   ofs_b<<"ASCII"<<std::endl;
   ofs_b<<"DATASET POLYDATA"<<std::endl;
   ofs_b<<"POINTS"<<" "<<points.size()*beads_ex.size()<<" "<<"float"<<std::endl;
   boost::format point("%f %f %f");
   for(size_t i=0,size=beads_ex.size();i<size;++i)
   {
      const Vector3D& center = beads_ex[i].first;
      for(size_t s=0,ssze=points.size();s<ssze;++s)
      {
         const Vector3D p = (points[s]*beads_ex[i].second+center);
         ofs_b<<point %(float)p.x %(float)p.y %(float)p.z<<std::endl;
      }
   }
   boost::format pole("3 %d %d %d");
   boost::format equator("4 %d %d %d %d");
   ofs_b<<"POLYGONS "<<(triangles.size()+quadrangles.size())*beads_ex.size()<<" "<<(4*triangles.size()+5*quadrangles.size())*beads_ex.size()<<std::endl;

   const int polygon_size = (points.size());
   for(size_t i=0,size=beads_ex.size();i<size;++i)
   {
      for(size_t t=0,stze=triangles.size();t<stze;++t)
      {
         ofs_b<<pole %(triangles[t].a+i*polygon_size) %(triangles[t].b+i*polygon_size) %(triangles[t].c+i*polygon_size)<<std::endl;
      }
      for(size_t q=0,sqze=quadrangles.size();q<sqze;++q)
      {
         ofs_b<<equator %(quadrangles[q].a+i*polygon_size) %(quadrangles[q].b+i*polygon_size) %(quadrangles[q].c+i*polygon_size) %(quadrangles[q].d+i*polygon_size)<<std::endl;
      }
   }

   if(!scalars.empty())
   {
      ofs_b<<"POINT_DATA "<<points.size()*beads_ex.size()<<std::endl;
      ofs_b<<"SCALARS "<<"beads "<<"float 1"<<std::endl;
      ofs_b<<"LOOKUP_TABLE "<<"default"<<std::endl;
      boost::format cl("%e ");
      for(size_t i=0, size=beads_ex.size(); i<size; ++i)
      {
         for(size_t j=0, sjze=points.size(); j<sjze; ++j)
         {
            ofs_b<<cl %(float)scalars[i];
         }
         ofs_b<<std::endl;
      }
   }
   ofs_b.close();
   return ;
}

void E2_base
(
   const unsigned int& step,
   const std::vector<Ellipsoid2D> E2s
)
{

   //tmp
   std::vector<Vector3D>   vertices;
   std::vector<Triangle>   triangles;
   std::vector<Quadrangle> quadrangles;

   const auto rot = [](const Vector3D& v,const double theta)
   {
      const Quaternion q(Vector3D(0.0,0.0,-1.0),(0.5*M_PI-theta));
      const auto res = q*v*bar(q);
      return Vector3D(res.b, res.c, res.d);
   };
   const auto apply = [&](const Ellipsoid2D& e2)->void
   {
      const double& x     = e2.x;
      const double& y     = e2.y;
      const double& theta = e2.theta;
      const double& a     = e2.a;
      const double& b     = e2.b;

      auto etmp = Ellipsoid(a,b,Devide_N);

      //rot
      std::vector<Vector3D>& ps = std::get<std::vector<Vector3D> >(etmp);
      std::for_each(ps.begin(),ps.end(),[&theta,&rot](auto& v){v=rot(v,theta);});
      //shift
      std::for_each(ps.begin(),ps.end(),[&x,&y](auto& v){v.x+=x;v.y+=y;});
      vertices    = ps;
      triangles   = std::get<std::vector<Triangle> > (etmp);
      quadrangles = std::get<std::vector<Quadrangle> > (etmp);
      return ;
   };
   try{
      apply(*E2s.begin());
   }catch(...){throw "empty";}

   const std::string str_step = boost::lexical_cast<std::string> (step);
   std::string e2d_file = output_file+"_e2d_"+str_step+".vtk";
   std::ofstream ofs_(e2d_file,std::ios::trunc);
   ofs_<<"# vtk DataFile Version 2.0"<<std::endl;
   ofs_<<"Ellipsoid 2D "<<input_file<<std::endl;
   ofs_<<"ASCII"<<std::endl;
   ofs_<<"DATASET POLYDATA"<<std::endl;
   ofs_<<"POINTS"<<" "<<E2s.size()*vertices.size()<<" "<<"float"<<std::endl;
   boost::format point("%f %f %f");
   for(int i=0,size=E2s.size();i<size;++i)
   {
      apply(E2s.at(i));
      for(int j=0,sjze=vertices.size();j<sjze;++j)
      {
         ofs_<<point %(float)vertices.at(j).x %(float)vertices.at(j).y %(float)vertices.at(j).z<<std::endl;
      }
   }
   boost::format pole("3 %d %d %d");
   boost::format equator("4 %d %d %d %d");
   ofs_<<"POLYGONS "<<(triangles.size()+quadrangles.size())*E2s.size()<<" "<<(4*triangles.size()+5*quadrangles.size())*E2s.size()<<std::endl;
   const int vsize = vertices.size();
   for(size_t i=0,size=E2s.size();i<size;++i)
   {
      for(size_t t=0,stze=triangles.size();t<stze;++t)
      {
         ofs_<<pole %(triangles[t].a+i*vsize) %(triangles[t].b+i*vsize) %(triangles[t].c+i*vsize)<<std::endl;
      }
      for(size_t q=0,sqze=quadrangles.size();q<sqze;++q)
      {
         ofs_<<equator %(quadrangles[q].a+i*vsize) %(quadrangles[q].b+i*vsize) %(quadrangles[q].c+i*vsize) %(quadrangles[q].d+i*vsize)<<std::endl;
      }
   }
}


void Beads
(
   const unsigned int& step,
   const std::vector<std::string>& list
)
{
   const std::vector<double> scalars = [&list]()
   {
      std::vector<double> res;
      if(!scalar_f){return res;}
      for(size_t i=0,i_size=list.size();i<i_size;++i)
      {
         std::vector<std::string> vs;
         boost::algorithm::split(vs, list.at(i), boost::is_any_of(" ,\t"));
      if(vs.size()==5)//index x y z scalar
      {
         res.push_back(boost::lexical_cast<double>(vs.back()));
      }
      if(vs.size()==7)//index qx qy qz px py pz
      {
         //scalar <- p^2
         res.push_back
            (
               std::pow(boost::lexical_cast<double>(vs.at(4)),2)+
               std::pow(boost::lexical_cast<double>(vs.at(5)),2)+
               std::pow(boost::lexical_cast<double>(vs.at(6)),2)
            );
      }
      else{std::cout<<"defect __LINE__"<<std::endl;exit(0);}
      }
      return res;
   }();
   return Beads_base(step,TextReader::cast_Coordinate(list),scalars);
}

void Beads_ex
(
   const unsigned int& step,
   const std::vector<std::string>& list
)
{
   const std::vector<double> scalars = [&list]()
   {
      std::vector<double> res;
      if(!scalar_f){return res;}
      for(size_t i=0,i_size=list.size();i<i_size;++i)
      {
         std::vector<std::string> vs;
         boost::algorithm::split(vs, list.at(i), boost::is_any_of(" ,\t"));
         if(vs.size()==6)//index x y z radius scalar
         {
            res.push_back(boost::lexical_cast<double>(vs.back()));
         }
         else{std::cout<<"defect __LINE__"<<std::endl;exit(0);}
      }
      return res;
   }();
   return Beads_ex_base(step,TextReader::cast_Beads_ex(list),scalars);
}


void Membrane_base
(
   const unsigned int& step,
   const unsigned int& shf_mode,
   const std::vector<Vector3D>& vertices,
   const std::vector<Triangle>& triangles,
   const std::vector<double>& scalars
)
{
   const std::string str_step = boost::lexical_cast<std::string> (step);
   const std::string str_mode = boost::lexical_cast<std::string> (shf_mode);
   std::string shf_file = output_file+"_shf_"+str_mode+"_"+str_step+".vtk";
   std::ofstream ofs_m(shf_file,std::ios::trunc);
   ofs_m<<"# vtk DataFile Version 2.0"<<std::endl;
   ofs_m<<"Mesh "<<input_file<<std::endl;
   ofs_m<<"ASCII"<<std::endl;
   ofs_m<<"DATASET POLYDATA"<<std::endl;
   ofs_m<<"POINTS"<<" "<<vertices.size()<<" "<<"float"<<std::endl;
   boost::format point("%e %e %e");
   for(int i=0,size=vertices.size();i<size;++i)
   {
      ofs_m<<point %(float)vertices[i].x %(float)vertices[i].y %(float)vertices[i].z<<std::endl;
   }
   ofs_m<<"POLYGONS"<<" "<<triangles.size()<<" "<<(4*triangles.size())<<std::endl;
   boost::format triangle("3 %d %d %d");
   for(size_t i=0,size=triangles.size();i<size;++i)
   {
      ofs_m<<triangle %triangles[i].a %triangles[i].b %triangles[i].c<<std::endl;
   }
   if(!scalars.empty())
   {
      ofs_m<<"POINT_DATA "<<vertices.size()<<std::endl;
      ofs_m<<"SCALARS "<<"beads "<<"float 1"<<std::endl;
      ofs_m<<"LOOKUP_TABLE "<<"default"<<std::endl;
      boost::format cl("%e ");
      for(size_t i=0, size=vertices.size(); i<size; ++i)
      {
         ofs_m<<cl %(float)scalars[i]<<std::endl;
      }
   }
   ofs_m.close();
}

void E2
(
   const unsigned int& step,
   const std::vector<std::string>& list
)
{
   return E2_base(step,TextReader::cast_Ellipsoid2D(list));
}
#endif
