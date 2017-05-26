#include "brillouin_zone.hpp"

brillouin_zone::brillouin_zone(const char* type)
{
  assert (std::strcmp(type,"square") == 0 && "wrong bz-type");
  Vector2d kk;
  kk << M_PI, 0.0;
	um_ = calc_poly(kk, 4);
  kk << M_PI, M_PI;
  bz_ = calc_poly(kk, 4);
}

brillouin_zone::brillouin_zone(const brillouin_zone &other)
{
	um_ = other.um_;
	bz_ = other.bz_;
}

brillouin_zone& brillouin_zone :: operator= (const brillouin_zone &other)
{
	// avoid self-assignment
	if ( this != &other )
  {       
		um_ = other.um_;
		bz_ = other.bz_;
  }
  return *this;
}

Polygon brillouin_zone::calc_poly(Vector2d kk, int num_rot)
{
  Polygon poly;
	Matrix2d rot;
  assert(num_rot == 2 || num_rot == 4 || num_rot == 6 && "invalid rot symmetry in calc_poly");
  double ang = 2*M_PI/num_rot;

	rot << cos(ang), -sin(ang),
	       sin(ang),  cos(ang);

	for ( int i = 0; i < num_rot; i = i + 1 )
	{
  	poly.push_back(Point(kk(0),kk(1)));
		kk = rot*kk;
	}
  return poly;
}

double brillouin_zone::dr2dx(double r, double ang) const
{
	Vector3d res_vec, dvec;
	Vector3d kvec(get_K_pt(ang));
	Vector3d uvec(get_U_pt(ang));
	dvec  = kvec - uvec;
	double unorm = uvec.norm();
	double dnorm = dvec.norm();
	double knorm = unorm + dnorm;
  assert( 0.0 <= r && r <= knorm && "invalid radius in dr2dx");

	if (0.0 <= r && r <= unorm)
  {
		return r;
  }
  else if (unorm <= r && r <= knorm)
  {
		return knorm - r;
  }
	else
  {
		return 0.0;
	}
}

Vector3d brillouin_zone::mapin(Vector3d kk) const
{	
	Vector3d fold_vec;
	double kx,ky;
	kx = kk(0);
	ky = kk(1);
	while (abs(kx) > M_PI)
	{
		if (kx < 0.0)
		{
			kx = kx + 2.0*M_PI;
		}
		else
		{
			kx = kx - 2.0*M_PI;
		}
	}

	while (abs(ky) > M_PI)
	{
		if (ky < 0.0)
		{
			ky = ky + 2.0*M_PI;
		}
		else
		{
			ky = ky - 2.0*M_PI;
		}
	}
	fold_vec << kx, ky, 0.0;
	return fold_vec; 
}

Vector3d brillouin_zone::r2xy(double r, double ang, Vector3d pt) const
{
	Vector3d res_vec, dvec;
	Vector3d kvec(get_K_pt(ang));
	Vector3d uvec(get_U_pt(ang));
	dvec  = kvec - uvec;
	double unorm = uvec.norm();
	double dnorm = dvec.norm();
	double knorm = unorm + dnorm;
  assert(0.0 <= r && r <= knorm && "invalid radius in r2xy");

	if (0.0 <= r && r <= unorm)
	{
		res_vec(0) = pt(0) + r*cos(ang);
		res_vec(1) = pt(1) + r*sin(ang);
		res_vec(2) = pt(2);
	}
	else if (unorm < r && r <= knorm)
	{
		res_vec(0) = pt(0) + uvec(0) + dvec(0)*(r - unorm)/dnorm;
		res_vec(1) = pt(1) + uvec(1) + dvec(1)*(r - unorm)/dnorm;
		res_vec(2) = pt(2);
	}
	else
	{
		res_vec(0) = 0.0;
		res_vec(1) = 0.0;
		res_vec(2) = 0.0;
	}
	return res_vec;
}

double brillouin_zone::get_patch_len(double phi) const
{
	Vector3d uvec(get_U_pt(phi));
  Vector3d dvec;
  dvec = get_K_pt(phi) - get_U_pt(phi);
	return uvec.norm() + dvec.norm();
}

Vector3d brillouin_zone::get_U_pt(double phi) const
{
	CGAL::Object result;
	Point pt;
  Direct  d1(cos(phi),sin(phi));
  Ray     r1(Point(0.0,0.0),d1);
  Vector3d  kres;

 	for (EdgeIterator ei = um_.edges_begin(); ei != um_.edges_end(); ++ei)
 	{
 		result = intersection(r1,*ei);
		if (CGAL::assign(pt,result)) 
 		{
   		kres(0) =  pt.x(); 
			kres(1) =  pt.y();
      kres(2) =  0.0;
 		}
 	}      
	return kres;
}

double brillouin_zone::angle_dist(double phi1, double phi2) const
{
  // return angle distance modulo pi
  while (abs(phi1) > M_PI) phi1 = phi1 - copysign(1.0,phi1)*2*M_PI;
  while (abs(phi2) > M_PI) phi2 = phi2 - copysign(1.0,phi2)*2*M_PI;

	return abs(phi1 - phi2 - copysign(1.0, phi1 - phi2)*2*M_PI*floor(abs(phi1 - phi2)/M_PI));
}

Vector3d brillouin_zone::get_K_pt(double phi) const
{
  double dist = 100;
  double dist_temp;
  Point pt;
  Vector3d kk(0.0,0.0,0.0);

 	for (VertexIterator vi = bz_.vertices_begin(); vi != bz_.vertices_end(); ++vi)
  {
		pt = *vi;
		dist_temp = angle_dist(atan2(pt.y(),pt.x()),phi);
    if (dist_temp < dist)
		{
			dist = dist_temp;
			kk << pt.x(), pt.y(), 0.0;
		}
	}
	return kk;
}

Polygon brillouin_zone::get_bz() const
{
	return bz_;
}
  
Polygon brillouin_zone::get_um() const
{
	return um_;
}

void brillouin_zone::print_discrt(const char* filename, int patch_num, int rad_num) const
{
   ofstream file_output(filename);
   int width1 = 18;
   double ang, rlim, r;
   Vector3d pt(0.0,0.0,0.0), kvec;
   file_output.precision(10); 
   file_output << setw(width1) << right << "angle" << setw(width1) << right << "rr";
   file_output << setw(width1) << right << "kx" << setw(width1) << right << "ky";  
   file_output << "\n";

   for (int nn = 0; nn < patch_num; nn++)
   {
		ang  = (2*nn + 1)*M_PI/patch_num; 
		rlim = get_patch_len(ang);
	  for (int nr = 0; nr < rad_num; nr++)
	  { 
			r = nr*rlim/rad_num;
			kvec = r2xy(r,ang,pt);
			file_output << setw(width1) << right << ang;
			file_output << setw(width1) << right << r;
			file_output << setw(width1) << right << kvec(0);
			file_output << setw(width1) << right << kvec(1);
      file_output << "\n";
		}
	}
}
