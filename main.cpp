bool scalar_f=false;
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <thread>
#include <cmath>
#include <queue>
#include "Vector3D.hpp"
#include "Quaternion.hpp"
#include "VTK_Geometry.hpp"
#include "TextReader.hpp"
#include "ReadFile.hpp"
#include "PBC.hpp"
#include "Sphere.hpp"
#include "Output.hpp"
#define RUN_NAME "YUBA"
#include "Logger.hpp"

std::pair<std::vector<Coordinate>,std::vector<Triangle> > Sphere_template;

class tag_value
{
   public:
   std::string tag;
   std::string value;
   tag_value(const std::string& v1, const std::string& v2){tag=v1;value=v2;}
};

class setting
{
   public:
      bool input_f;
      bool output_f;
      bool step_f;
      bool mode_f;
      setting()
      {
         input_f=false; output_f=false; step_f=true; mode_f=false;
      }
      bool is_fail() const
      {
         return (input_f==false||output_f==false||step_f==false||mode_f==false)?true:false;
      }
};

class active_step_unit
{
   public:
      unsigned int bgn;
      unsigned int end;
      active_step_unit(){ bgn=INT_MAX, end=INT_MAX; }
      active_step_unit(const unsigned int& a)
      {
         bgn=a;
         end=a;
      }
      active_step_unit(const unsigned int& a, const unsigned int& b)
      {
         bgn=a;
         end=b;
      }
      bool on_range(const unsigned int& val) const
      {
         return ( (bgn<=val) && (val<=end) ) ? true : false ;
      }
};

std::vector<std::string>str_list(Getline& stream, const std::string& tag, const std::string& str_step);

/////////////////
std::vector<tag_value>            tags;
std::vector<active_step_unit>    steps;

std::string                       mode;
std::set<unsigned int>   steps_in_file;

Getline                vertices_stream;
Getline               triangles_stream;
Getline                   beads_stream;
Getline                   bonds_stream;
Getline             bond_points_stream;
Getline                      E2_stream;
Getline                      E3_stream;

double                          Tube_R;
int                            thd =-1;
const Vector3D             null_vector;
std::vector<unsigned int>  shf_ns;

/////////////////
void step_picker (const std::string& req,bool& f);
void Membrane
(
   const unsigned int& step,
   const std::vector<std::string>& list_vertices,
   const std::vector<std::string>& list_triangles
);
void Bonds
(
   const unsigned int& step,
   const std::vector<std::string>& list_beads,
   const std::vector<std::string>& list_bonds
);

void help();
void help(const std::string& additional, std::string message="none");

int main(int argc, char* argv[])
{
   //argv -> argments vector
   std::vector<std::string> argments;
   setting flags;
   for(int i=0;i<argc;++i){ argments.push_back(argv[i]); }

   //argments ---split and cast---> tags
   for(size_t i=0,size=argments.size();i<size;++i)
   {
      std::vector<std::string> vsplit;
      boost::algorithm::split(vsplit,argments[i],boost::is_any_of("="));
      if(vsplit.size()==2)
      {
         tags.push_back(tag_value(vsplit[0],vsplit[1]));
      }
      //help
      if("-h"==argments[i]){help();exit(1);}
   }

   //setting
   for(size_t i=0,size=tags.size();i<size;++i)
   {
      if(tags[i].tag=="input" ||
         tags[i].tag=="in"      ){ input_file=tags[i].value;   flags.input_f =true; }
      if(tags[i].tag=="output"||
         tags[i].tag=="out"     ){ output_file=tags[i].value;  flags.output_f=true; }
      if(tags[i].tag=="mode")
      {
         const std::string& md = tags[i].value;
         if
         (
            ("volvox"==md)||("tetra"==md)||("membrane"==md)||("beads"==md)||("bonds"==md)||("beads_ex"==md)||
            ("ellipsoid2d"==md)||("ellipsoid3d==md")
         )
         {
            mode=md;
            flags.mode_f=true;
         }
      }
      try{
         
         if(boost::iequals(tags[i].tag, "Bead_SIZE"))
         { 
            Bead_SIZE=boost::lexical_cast<double>(tags[i].value);
         }
         
         if(boost::iequals(tags[i].tag, "Bond_SIZE"))
         { Tube_R=boost::lexical_cast<double>(tags[i].value);}

         if(boost::iequals(tags[i].tag, "thread") || boost::iequals(tags[i].tag, "th"))
         { thd=boost::lexical_cast<int> (tags[i].value);}

         if(boost::iequals(tags[i].tag, "Bead_Grain") || boost::iequals(tags[i].tag, "bg") || 
            boost::iequals(tags[i].tag,"Grain"))
         {
            if(boost::iequals(tags[i].value,"finest")){ Devide_N= 4;}
            if(boost::iequals(tags[i].value,"fine")){ Devide_N= 3;}
            if(boost::iequals(tags[i].value,"normal")){ Devide_N= 2;}
            if(boost::iequals(tags[i].value,"rough")){ Devide_N= 1;}
            if(boost::iequals(tags[i].value,"roughest")){ Devide_N= 0;}
         }
         if(boost::iequals(tags[i].tag, "Bond_Grain") || boost::iequals(tags[i].tag, "bog"))
         {
            if(boost::iequals(tags[i].value,"finest")){ Devide_tube= 100;}
            if(boost::iequals(tags[i].value,"fine")){ Devide_tube= 70;}
            if(boost::iequals(tags[i].value,"normal")){ Devide_tube= 50;}
            if(boost::iequals(tags[i].value,"rough")){ Devide_tube= 30;}
            if(boost::iequals(tags[i].value,"roughest")){ Devide_tube= 10;}
         }
         if(boost::iequals(tags[i].tag, "scalar") || boost::iequals(tags[i].tag, "sc"))
         {
            if(boost::iequals(tags[i].value,"true"))
            {
               scalar_f=true;
            }
         }
         if(boost::iequals(tags[i].tag,"PBC"))
         {
            std::vector<std::string> vsplit_pbc;
            boost::algorithm::split(vsplit_pbc,tags[i].value,boost::is_any_of(","));
            for(size_t i_=0; i_<vsplit_pbc.size(); ++i_)
            {
               if(boost::iequals(vsplit_pbc[i_],"x")){ PBC_f[0]=true; }
               if(boost::iequals(vsplit_pbc[i_],"y")){ PBC_f[1]=true; }
               if(boost::iequals(vsplit_pbc[i_],"z")){ PBC_f[2]=true; }
            }
         }
      }
      catch(boost::bad_lexical_cast &)
      {
         help("Arguments fail");
      }
   }
   //default
   if(Bead_SIZE<=0){Bead_SIZE=1.0;}
   if(Tube_R<=0){Tube_R=0.5;}

   Sphere_template=Sphere(0.5*(boost::lexical_cast<double> (Bead_SIZE)));
   Sphere_template_ex=Regular_Sphere(Devide_N*20);

   //extract steps
   std::ifstream ifs(input_file);
   if(ifs.fail()){help("Input file fail",input_file+" not exist.");}
   std::string reading;
   bool cell_collect_f=false;
   while(getline(ifs,reading))
   {
      std::vector<std::string> vs;
      boost::algorithm::split(vs,reading,boost::is_any_of(" ,\t"));
      if(vs[0]=="#")
      {
         if(vs[1]!="Cell")//normal sequence
         {
            steps_in_file.insert(boost::lexical_cast<unsigned int> (vs[2]));
         }
         else if(vs[1]=="Cell")
         {  
            if(vs.size()>3)
            {
               Cell_length[0]=boost::lexical_cast<double> (vs[2]);
               Cell_length[1]=boost::lexical_cast<double> (vs[3]);
               Cell_length[2]=boost::lexical_cast<double> (vs[4]);
            }
            else
            {
               Cell_length[0]=boost::lexical_cast<double> (vs[2]);
               Cell_length[1]=boost::lexical_cast<double> (vs[2]);
               Cell_length[2]=boost::lexical_cast<double> (vs[2]);
            }
            cell_collect_f=true;
         }
      }
      else if(cell_collect_f)
      {
         //string -> vector of center
         const Vector3D center(
            boost::lexical_cast<double> (vs[0]),
            boost::lexical_cast<double> (vs[1]),
            boost::lexical_cast<double> (vs[2])
         );
         center_of_cell=center;
         cell_collect_f=false;
         cell_f=true;
      }
   }
   ifs.close();


   for(int i=0,size=tags.size();i<size;++i)
   {  //setting: step regions
      if(boost::iequals(tags[i].tag, "step"))
      {
         flags.step_f=false;
         step_picker(tags[i].value,flags.step_f);
      }
   }

   if(steps_in_file.empty()){flags.step_f=false;}
   if(flags.is_fail()){help("Essential argment(s)");}

   vertices_stream.set(input_file);
   triangles_stream.set(input_file);
   beads_stream.set(input_file);
   bonds_stream.set(input_file);
   bond_points_stream.set(input_file);
   E2_stream.set(input_file);
   E3_stream.set(input_file);

   if(cell_f)
   {  //Cell for the boundary conditions
      std::vector<Vector3D> points;
      const Vector3D& c = center_of_cell;
      for(int x=0;x<=1;++x)
      {
         for(int y=0;y<=1;++y)
         {
            for(int z=0;z<=1;++z)
            {
               points.push_back
               (
                  Vector3D (
                     c.x + ((x==0)?Cell_length[0]/2.0:-Cell_length[0]/2.0),   
                     c.y + ((y==0)?Cell_length[1]/2.0:-Cell_length[1]/2.0),   
                     c.z + ((z==0)?Cell_length[2]/2.0:-Cell_length[2]/2.0)
                  )
               );
            }
         }
      }
      std::string cell_file = output_file+"_cell.vtk";
      std::ofstream ofs_c(cell_file,std::ios::trunc);
      ofs_c<<"# vtk DataFile Version 2.0"<<std::endl;
      ofs_c<<"CELL "<<input_file<<std::endl;
      ofs_c<<"ASCII"<<std::endl;
      ofs_c<<"DATASET POLYDATA"<<std::endl;
      ofs_c<<"POINTS"<<" "<<points.size()<<" "<<"float"<<std::endl;
      boost::format point("%f %f %f");

      //Faces of the cube
      for(int i=0,size=points.size();i<size;++i)
      {
         ofs_c<<point %points[i].x %points[i].y %points[i].z<<std::endl;
      }
      ofs_c<<"POLYGONS"<<" "<<"6"<<" "<<"30"<<std::endl;
      ofs_c<<"4 0 1 5 4"<<std::endl;
      ofs_c<<"4 0 2 3 1"<<std::endl;
      ofs_c<<"4 5 1 3 7"<<std::endl;
      ofs_c<<"4 4 6 7 5"<<std::endl;
      ofs_c<<"4 0 2 6 4"<<std::endl;
      ofs_c<<"4 6 2 3 7"<<std::endl;
      ofs_c.close();
   }

   //selecting of target steps
   std::vector<unsigned int> request;
   for(auto i=steps_in_file.begin(); i!=steps_in_file.end(); ++i)
   {
      if(!steps.empty())
      {
         for(int s=0,ssze=steps.size();s<ssze;++s)
         {
            if(steps[s].on_range(*i))
            {
               request.push_back(*i);
            }
         }
      }
      else
      {
         request.push_back(*i);
      }
   }

   try
   {
      #define JOB_N_DEF 2
      const int thd_max=std::thread::hardware_concurrency()-2;
      if(thd==-1){thd=thd_max;}
      const int JOB_N=std::max(1,std::min(thd,thd_max));
      for(int i=0,size=request.size();i<size;i+=JOB_N)
      {
         const int iJ = i+JOB_N;
         const int& thlimit = (iJ>size)?size:iJ;
         if(boost::iequals("tetra", mode)||boost::iequals("bonds", mode))
         {//WAIDE
            std::vector<std::thread>                       threads;
            std::vector<std::vector<std::string> >  str_list_beads;
            std::vector<std::vector<std::string> >  str_list_bond_points;
            std::vector<std::vector<std::string> >  str_list_bonds;
            std::string Bond_Point_type="none";
            for(int j=0,sjze=tags.size();j<sjze;++j)
            {
               if(boost::iequals(tags[j].tag, "Bond_Point_type") || boost::iequals(tags[j].tag, "BPT"))
               {
                  Bond_Point_type=tags[j].value;
                  break;
               }
            }
            if(Bond_Point_type=="none"){Bond_Point_type="Coordinate";}
            
            for(int t=i;t<thlimit;++t)
            {
               str_list_bond_points.push_back
               (
                  str_list(bond_points_stream, Bond_Point_type, boost::lexical_cast<std::string>(request[t]))
               );
               str_list_bonds.push_back
               (
                  str_list(bonds_stream, "Bond", boost::lexical_cast<std::string>(request[t]))
               );
            }

            for(int t=i;t<thlimit;++t)
            {
               threads.push_back(std::thread([&,t]()
               {
                  Bonds(request[t], str_list_bond_points[t-i],str_list_bonds[t-i]);
               }));
            }

            std::for_each(threads.begin(),threads.end(),[](std::thread& t){t.join();});
         }
         if(boost::iequals("membrane",mode)||boost::iequals("volvox",mode))
         {//HIBIKI
            std::vector<std::thread>                          threads;
            std::vector<std::vector<std::string> >  str_list_vertices;
            std::vector<std::vector<std::string> > str_list_triangles;
            std::vector<Vector3D>                             centers;
            for(int t=i;t<thlimit;++t)
            {
               str_list_vertices.push_back
               (
                  str_list(vertices_stream,"Coordinate",boost::lexical_cast<std::string>(request[t]))
               );
               str_list_triangles.push_back
               (
                  str_list(triangles_stream,"Triangle",boost::lexical_cast<std::string>(request[t]))
               ); 
            }
            for(int t=i;t<thlimit;++t)
            {
               threads.push_back(std::thread([&,t]()
               {
                  Membrane(request[t],str_list_vertices[t-i],str_list_triangles[t-i]);
               }));
            }
            std::for_each(threads.begin(),threads.end(),[](std::thread& t){t.join();});
         }
         if(boost::iequals("beads",mode)||boost::iequals("volvox",mode)||boost::iequals("tetra", mode))
         {//BOTH
            std::vector<std::thread>                       threads;
            std::vector<std::vector<std::string> >  str_list_beads;
            for(int t=i;t<thlimit;++t)
            {
               str_list_beads.push_back(str_list(beads_stream,"Bead",boost::lexical_cast<std::string>(request[t])));
            }
            for(int t=i;t<thlimit;++t)
            {
               threads.push_back(std::thread([&,t]()
               {
                  Beads(boost::lexical_cast<unsigned int>(request[t]),str_list_beads[t-i],scalar_f);
               }));
            }
            std::for_each(threads.begin(),threads.end(),[](std::thread& t){t.join();});
         
         }
         if(boost::iequals("beads_ex",mode))
         {
            std::vector<std::thread>                 threads;
            std::vector<std::vector<std::string> >   str_list_beads_ex;
            for(int t=i;t<thlimit;++t)
            {
               str_list_beads_ex.push_back(str_list(beads_stream,"Bead_ex",boost::lexical_cast<std::string>(request[t])));
            }
            for(int t=i;t<thlimit;++t)
            {
               threads.push_back
                  (
                     std::thread
                     (
                        [&,t]() 
                        {
                           Beads_ex(boost::lexical_cast<unsigned int>(request[t]),str_list_beads_ex[t-i],scalar_f);
                        }
                     ) 
                  );
            }
            std::for_each(threads.begin(),threads.end(),[](std::thread& t){t.join();});
         }
         if(boost::iequals("ellipsoid2d",mode))
         {
            std::vector<std::thread>                 threads;
            std::vector<std::vector<std::string> >   str_list_e2;
            for(int t=i;t<thlimit;++t)
            {
               str_list_e2.push_back(str_list(E2_stream,"Ellipsoid2D",boost::lexical_cast<std::string>(request[t])));
            }
            for(int t=i;t<thlimit;++t)
            {
               threads.push_back
                  (
                     std::thread
                     (
                        [&,t]() 
                        {
                           E2(boost::lexical_cast<unsigned int>(request[t]),str_list_e2[t-i]);
                        }
                     ) 
                  );
            }
            std::for_each(threads.begin(),threads.end(),[](std::thread& t){t.join();});
            
         }
      }
   }
   catch(const std::exception& ex)
   {
      std::cout<<ex.what()<<std::endl;
      help("FATAL ERROR: Please mail me!");
   }

   return EXIT_SUCCESS;
}

void Membrane
(
   const unsigned int& step,
   const std::vector<std::string>& list_vertices,
   const std::vector<std::string>& list_triangles//,
)
{
   const std::string str_step = boost::lexical_cast<std::string> (step);

   //READ and CAST info. of "step"
   std::vector<Coordinate> vertices       = TextReader::cast_Coordinate(list_vertices);
   const std::vector<Triangle>& triangles = TextReader::cast_Triangle(list_triangles);

   //Mesh Output
   std::string mesh_file = output_file+"_mesh_"+str_step+".vtk";
   std::ofstream ofs_m(mesh_file,std::ios::trunc);
   ofs_m<<"# vtk DataFile Version 2.0"<<std::endl;
   ofs_m<<"Mesh "<<input_file<<std::endl;
   ofs_m<<"ASCII"<<std::endl;
   ofs_m<<"DATASET POLYDATA"<<std::endl;
   ofs_m<<"POINTS"<<" "<<vertices.size()<<" "<<"float"<<std::endl;
   boost::format point("%e %e %e");
   for(size_t i=0,size=vertices.size();i<size;++i)
   {
      ofs_m<<point %(float)vertices[i].v.x %(float)vertices[i].v.y %(float)vertices[i].v.z<<std::endl;
   }
   ofs_m<<"POLYGONS"<<" "<<triangles.size()<<" "<<(4*triangles.size())<<std::endl;
   boost::format triangle("3 %d %d %d");
   for(size_t i=0,size=triangles.size();i<size;++i)
   {
      ofs_m<<triangle %triangles[i].a %triangles[i].b %triangles[i].c<<std::endl;
   }
   ofs_m.close();
}


void Bonds
(
   const unsigned int& step, 
   const std::vector<std::string>& list_vertices, 
   const std::vector<std::string>& list_bonds
)
{
   const std::string str_step = boost::lexical_cast<std::string> (step);
   const std::vector<Coordinate>& bonds_points = TextReader::cast_Coordinate(list_vertices);
   const std::vector<Bond>& bonds = TextReader::cast_Bond(list_bonds);
   
   typedef int scalar;
   std::vector<std::pair<std::vector<Vector3D>, scalar> > rectangles_info;

   for(int i=0, size=bonds.size(); i<size; ++i)
   {  
      std::vector<std::pair<Vector3D, Vector3D> > List_of_MakeBond;
      bool pending =true;

      Vector3D A=bonds_points[bonds[i].v1].v;
      const Vector3D B=bonds_points[bonds[i].v2].v;
      do
      {
         bool repeat_f=true;
         const Vector3D B_dash=Nearest_atom(A,B);
         std::vector<Vector3D> cross_point_list;
         for(int d=0, Dimention=3; d<Dimention; ++d)
         {
            if(PBC_f[d] && std::abs(A[d]-B[d])>0.5*Cell_length[d])
            {
               for(int j=0; j<2; ++j)
               {
                  Vector3D normal_vector_of_plane;
                  normal_vector_of_plane[d]=1.0;
                  const Vector3D result=[&](const Vector3D& a, const Vector3D& b)->Vector3D
                  {
                     const double parameter=(Cell_length[d]*j-(normal_vector_of_plane*a))
                                             /(normal_vector_of_plane*(b-a));
                     const Vector3D cross_point=(a+(b-a)*parameter);
                     return cross_point; 
                  }(A, B_dash);
                  cross_point_list.push_back(result);
               }
               repeat_f=false;   
            }
            else if(d==(Dimention-1)&& repeat_f)
            {
               List_of_MakeBond.push_back(std::pair<Vector3D, Vector3D>(A,B));
               pending=false;
            }
         }
         const Vector3D cross_point=[&]()
         {
            Vector3D result;
            double min_norm=DBL_MAX;
            for(int j=0, sjze=cross_point_list.size(); j<sjze; ++j)
            {
               if(!isnan(cross_point_list[j].x))
               {
                  if(min_norm>(A-cross_point_list[j]).norm())
                  {
                     result=cross_point_list[j];
                     min_norm=(A-cross_point_list[j]).norm();
                  }
               }
            }
            return result; 
         }();

         if(!repeat_f)
         {
            List_of_MakeBond.push_back(std::pair<Vector3D, Vector3D>(A, cross_point));
            const Vector3D cross_point_dash=Nearest_atom(B, cross_point);
            if(is_Box_in(cross_point_dash))
            {
               List_of_MakeBond.push_back(std::pair<Vector3D, Vector3D>(cross_point_dash, B));
               pending=false;
            }else
            {
               const Vector3D cross_point_dash2=[&](const Vector3D& P)
               {
                  Vector3D result;
                  Vector3D diff;
                  for(int d=0; d<3; ++d)
                  {
                     if(cross_point[d]==0.0 || cross_point[d]==Cell_length[d])
                     {
                        diff[d]+=Cell_length[d];
                     }
                  }
                  const Vector3D Plus  = cross_point+diff;
                  const Vector3D Minus = cross_point-diff;
                  if((Plus-P).norm()<(Minus-P).norm())
                  {   return Plus;  }
                  else
                  {   return Minus; }
               }(B);
               A=cross_point_dash2;
            }
         }
      }while(pending);

      for(int j=0, sjze=List_of_MakeBond.size(); j<sjze; ++j)
      {
         rectangles_info.push_back
         (
            std::pair<std::vector<Vector3D>, int>
            (make_Tube(List_of_MakeBond[j].first, List_of_MakeBond[j].second, Tube_R), bonds[i].scalar)
         );
      }
   }// bond for fin.

   std::vector<Rectangle> rectangles_index(Devide_tube);
   for(int i=0; i<Devide_tube; ++i)
   {
      rectangles_index[i].a = i;                 
      rectangles_index[i].b = (i+1)%Devide_tube; 
      rectangles_index[i].c = rectangles_index[i].a + Devide_tube;
      rectangles_index[i].d = rectangles_index[i].b + Devide_tube;
   }
   
   //Output 
   std::string tube_file = output_file+"_tube_"+ str_step+".vtk";
   std::ofstream ofs_t(tube_file, std::ios::trunc);
   ofs_t<<"# vtk DataFile Version 2.0"<<std::endl;
   ofs_t<<"Tube "<<std::endl;
   ofs_t<<"ASCII"<<std::endl;
   ofs_t<<"DATASET POLYDATA"<<std::endl;

   ofs_t<<"POINTS"<<" "<<rectangles_info.size()*rectangles_info[0].first.size()<<" "<<"float"<<std::endl;
   boost::format pt("%f %f %f");
   for(size_t i=0; i<rectangles_info.size(); ++i)
   {
      auto output=rectangles_info[i].first;
      for(size_t j=0, sjze=output.size(); j<sjze; ++j)
      {
         ofs_t<<pt %(float)output[j].x %(float)output[j].y %(float)output[j].z <<std::endl;
      }
   }

   boost::format pg("4 %d %d %d %d");
   ofs_t<<"POLYGONS"<<" "<<rectangles_index.size()*rectangles_info.size()
        <<" "<<(5*rectangles_index.size()*rectangles_info.size())<<std::endl;
   for(size_t i=0; i<rectangles_info.size(); ++i)
   {
      for(size_t t=0; t<rectangles_index.size(); ++t)
      {
         const int RectanglesSize_PerOneBond=rectangles_info[i].first.size();
         ofs_t<<pg   %(rectangles_index[t].a+i*RectanglesSize_PerOneBond)
                     %(rectangles_index[t].b+i*RectanglesSize_PerOneBond)
                     %(rectangles_index[t].d+i*RectanglesSize_PerOneBond)
                     %(rectangles_index[t].c+i*RectanglesSize_PerOneBond)<<std::endl;
      }
   }
   if(scalar_f)
   {
      ofs_t<<"CELL_DATA "<<rectangles_index.size()*rectangles_info.size()<<std::endl;
      ofs_t<<"SCALARS "<<"bonds "<<"float 1"<<std::endl;
      ofs_t<<"LOOKUP_TABLE "<<"default"<<std::endl;
      boost::format cl("%d ");
      for(size_t i=0, size=rectangles_info.size(); i<size; ++i)
      {
         for(size_t j=0; j<rectangles_index.size(); ++j)
         {
            ofs_t<<cl %rectangles_info[i].second<<std::endl;
         }
      }
   }
   ofs_t.close();
}

inline unsigned int str2step(const std::string& str)
{
   unsigned int result=INT_MAX;
   try{
      if(boost::iequals(str,"last"))
      {
         unsigned int tmp_max = 0;
         for(auto m=steps_in_file.begin();m!=steps_in_file.end();++m)
         {
            const unsigned int& tmp = *m;
            if(tmp_max<tmp){ tmp_max=tmp; }
         }
         result = tmp_max;
      }
      else
      {
         result = boost::lexical_cast<unsigned int> (str);
      }
   }
   catch(...)
   {
      help("unexpected step");
   }
   return result;
}


inline void step_picker (const std::string& req, bool& setuped_f)
{
   //split by comma
   std::vector<std::string> fst;
   const std::string fst_delimiter=",";
   boost::split(fst,req,boost::is_any_of(fst_delimiter));
   if(!fst.empty()){setuped_f=true;}

   if((4==fst.size())&&("..."==fst[2]))
   {//for Arithmetic progression
      const unsigned int first_step  = str2step(fst[0]);
      const unsigned int second_step = str2step(fst[1]);
      const unsigned int final_step  = str2step(fst[3]);
      if((second_step<first_step)&&(first_step<final_step)){help("Arithmetic progression failure");}
      if((first_step==second_step)){help("Arithmetic progression failure");}
      const long diff = second_step-first_step;
      long int   s_   = first_step;
      while(((diff>0)&&s_<=(long int)final_step)||((diff<0)&&(s_>=(long int)final_step)))
      {
         steps.push_back(active_step_unit(s_));
         s_+=diff;
      }
   }
   else
   {
      //split by - , and Each begin and end of the region is stored in snd[].
      for(int i=0,size=fst.size();i<size;++i)
      {
         std::vector<std::string> snd;
         const std::string snd_delimiter="-";
         boost::split(snd,fst[i],boost::is_any_of(snd_delimiter));
         
         if(snd.size()==1)
         {
            steps.push_back(active_step_unit(str2step(snd[0])));
         }
         else if(snd.size()==2)
         {
            steps.push_back(active_step_unit(str2step(snd[0]), str2step(snd[1])));
         }
      }
   }
}

std::vector<std::string>
str_list
(Getline& stream, const std::string& tag, const std::string& str_step)
{
   std::vector<std::string> results;
   bool read_f = false;
   bool fin_f  = false;

   while(stream.is_open()&&!fin_f)
   {
      try{
      const std::string tmp = stream.get();  
      std::vector<std::string> vs;
      boost::algorithm::split(vs,tmp,boost::is_any_of(" ,\t"));
   
      if("#"==vs[0])
      {
         if(read_f)
         {
            fin_f=true;
            stream.back();
         }
         else if((tag==vs[1])&&(str_step==vs[2]))
         {
            read_f=true;
         }
      }
      else
      {
         if(read_f)
         {
            if((int)vs.size()>1)
            {
               results.push_back(tmp);
            }
            else{fin_f=true;}
         }
      }
      }catch(...){}
   }
   if(results.empty())
   {
      std::cout<<tag<<" "<<str_step<<std::endl;
      std::cout<<"Error : str_list empty."<<std::endl;
      exit(1);
   }
   return results;
}

inline void help()
{
   printf("\x1b[35m");
   std::cout<<"                       Y U B A"<<std::endl;
   std::cout<<std::endl;
   printf("\x1b[32m");
   std::cout<<"                    Essential:"<<std::endl;
   std::cout<<"       Input File : input=huga.dat or in=huga.dat"<<std::endl;
   std::cout<<"      Output File : output=huga    or out=huga"<<std::endl;
   std::cout<<"        Data Mode : mode=piyo, ex. tetra, volvox, membrane, beads, beads_ex, bonds,"<<std::endl;
   std::cout<<"                                   ellipsoid2d, ellipsoid3d"<<std::endl;
   std::cout<<"   If you select bonds or tetra mode, add the option 'Bond_Point_type'. "<<std::endl;
   std::cout<<std::endl;
   printf("\x1b[36m");
   std::cout<<"                       Optional:"<<std::endl;
   std::cout<<"          Steps   : step=a, step=X-Y,A-B,W,L , and step=a,b,...,n ( default: all steps )"<<std::endl;
   std::cout<<" Num. of threads  : thread=6 or th=6 ( default: 2 )"<<std::endl;
   std::cout<<"   Beads Diamter  : Bead_SIZE=1.23   ( default: 1.0 )"<<std::endl;
   std::cout<<"   Bonds Diamter  : Bond_SIZE=1.23   ( default: 0.5 )"<<std::endl;
   std::cout<<"     Bond Scalar  : scalar=true or sc=true ( default: false )"<<std::endl;
   std::cout<<"             PBC  : PBC=x PBC=x,y,z  ( default: Free Boundary Condition ) "<<std::endl;
   std::cout<<" Beads Roughness  : Bead_Grain=ho or BG=ho  [rough, normal, fine]"<<std::endl;
   std::cout<<" Bonds Roughness  : Bond_Grain=ho or BoG=ho [rough, normal, fine]"<<std::endl;
   std::cout<<"Bonds Point Type  : Bond_Point_type=type (default: Coordinate)"<<std::endl;
   std::cout<<std::endl;
   printf("\x1b[31m");
   printf("\x1b[39m"); 
   printf("\x1b[49m");
   std::cout<<"       Bug report : @master-yde ( tetra, bond, scalar )"<<std::endl;
   std::cout<<"                    @misteltein ( otherwise )"<<std::endl;
}
inline void help(const std::string& additional, std::string message)
{
   help();
   std::cout<<"\nComments: "<<additional<<std::endl;
   std::cout<<"Message: "<<message<<std::endl;
   exit(1);
}
