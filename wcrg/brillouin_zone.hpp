#ifndef BRILLOUIN_ZONE
#define BRILLOUIN_ZONE

#include "general_header.hpp"

class brillouin_zone
{
  private:
		Polygon um_;
    Polygon bz_;

  public:
    brillouin_zone(const char* type);
    brillouin_zone(const brillouin_zone &other);
	  brillouin_zone& operator=(const brillouin_zone &other);
 
		double dr2dx(double r, double ang) const; //integration measure
		Vector3d r2xy(double r, double ang, Vector3d pt) const; //radial to kartesian coordinates
		Vector3d mapin(Vector3d kk) const; //fold back to BZ

  	Vector3d get_U_pt(double phi) const;
  	Vector3d get_K_pt(double phi) const;
  	double get_patch_len(double phi) const;
  	double angle_dist(double phi1, double phi2) const;
  	Polygon calc_poly(Vector2d kk, int num_rot);

    Polygon get_bz() const;
    Polygon get_um() const;

		void print_discrt(const char* filename, int patch_num, int rad_num) const; //print entire discretiation paths
};

#endif
