/*    Copyright (C) 2012 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */
#ifndef fslsurface_STRUCTS_H
#define fslsurface_STRUCTS_H

#ifdef __linux
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glext.h>
#else
#include <OpenGL/gl.h>
#endif
#include <cmath>

namespace fslsurface_name{
	enum MarchingCubesMode {EQUAL_TO, GREATER_THAN, GREATER_THAN_EQUAL_TO};

    struct float4{
        float4(float xin, float yin, float zin, float win)
        {
            x=xin;
            y=yin;
            z=zin;
            w=win;

        }
        float4()
        {

        }
        float x,y,z,w;
    };

	template<class T>
	struct vec3 {
        vec3(const vec3<T> & v)
        {
            x=v.x;
            y=v.y;
            z=v.z;
        }
		vec3(T xin, T yin, T zin)
        {
            x=xin;
            y=yin;
            z=zin;
        }
        vec3()
        {
			
        }
        vec3<T>& operator-=( const vec3<T> &other );
        vec3<T>& operator+=( const vec3<T> &other );

        vec3<T>& operator*=( const vec3<T> &other );
        vec3<T>& operator*=( const T & sc );
        vec3<T>& operator/=( const T & sc );

        const vec3<T> operator*(const T &other) const;
        const vec3<T> operator+(const vec3<T> &other) const;
        const vec3<T> operator-(const vec3<T> &other) const;

        void normalize();
        double norm();
        T x,y,z;
	};
    
    template<class T>
    void vec3<T>::normalize()
    {
        T l2 = (this->x * this->x) +  (this->y * this->y) + (this->z * this->z) ;
        
        *this /= sqrt(l2);
    }
    
    template<class T>
    double vec3<T>::norm()
    {
        T l2 = (this->x * this->x) +  (this->y * this->y) + (this->z * this->z) ;
        return static_cast<double>(sqrt(l2));
    }
    
    template<class T>
    vec3<T>& vec3<T>::operator-=( const vec3<T> & v1 )
    {
        this->x -= v1.x;
        this->y -= v1.y;
        this->z -= v1.z;

        return *this;
    }
    
    template<class T>
    vec3<T>& vec3<T>::operator+=( const vec3<T> & v1 )
    {
        this->x += v1.x;
        this->y += v1.y;
        this->z += v1.z;
        
        return *this;
    }
    
    template<class T>
    vec3<T>& vec3<T>::operator*=( const vec3<T> & v1 )
    {
        this->x *= v1.x;
        this->y *= v1.y;
        this->z *= v1.z;
        
        return *this;
    }
    template<class T>
    vec3<T>& vec3<T>::operator*=( const T & sc )
    {
        this->x *= sc;
        this->y *= sc;
        this->z *= sc;
        
        return *this;
    }
    template<class T>
    vec3<T>& vec3<T>::operator/=( const T & sc )
    {
        this->x /= sc;
        this->y /= sc;
        this->z /= sc;
        
        return *this;
    }
    template<class T>
    const vec3<T> vec3<T>::operator+(const vec3<T> &other) const {
        return vec3<T>(*this) += other;           
    }
    template<class T>
    const vec3<T> vec3<T>::operator-(const vec3<T> &other) const {
        return vec3<T>(*this) -= other;           
    }
    template<class T>
    const vec3<T> vec3<T>::operator*(const T &sc) const {
        return vec3<T>(*this) *= sc;           
    }

  /*  
    template<class T>
    vec3<T> subtract( const vec3<T> & v1, const vec3<T> & v2 )
    {
        return vec3<T>(v1.x-v2.xm,v1.y-v2.y,v1.z - v2.z);
    }
    */
	
	template <class T>
	struct vertex {
		
		vertex(T xin, T yin, T zin){
			x=xin;
			y=yin;
			z=zin;
			nx = 0;
			ny = 0;
			nz = 0;
			sc=0;
		}	
		vertex(T sc_in){
			x = 0;
			y = 0;
			z = 0;
			nx = 0;
			ny = 0;
			nz = 0;
			sc=sc_in;
		}
		vertex(){
			
		}
		
        vec3<T> coord() const;
      
		T x;
		T y;
		T z;
		T nx;
		T ny;
		T nz;
		T sc;
	};

    template<class T>
    vec3<T> vertex<T>::coord() const
    {
        return vec3<T>(this->x,this->y,this->z);
    }
    
    
    struct float3{
        float3(float xin, float yin, float zin)
        {
            x=xin;
            y=yin;
            z=zin;
        }
        float3()
        {

        }
        float x,y,z;
    };

    
    struct float2{
        float2(float xin, float yin)
        {
            x=xin;
            y=yin;
        }
        float2()
        {

        }
        float x,y;
    };
    struct bool3{
        bool3(bool xin, bool yin, bool zin)
        {
            x=xin;
            y=yin;
            z=zin;
        }
        bool3()
        {

        }
        bool x,y,z;
    };



    struct int2{
        int2(int xin,int yin)
        {
            x=xin;
            y=yin;
        }
        int2()
        {

        }

        float x,y;
    };

    struct material{

        material(float  ar,float  ag,float  ab,float  aa, \
                 float  dr,float  dg,float  db,float  da, \
                 float  sr,float  sg,float  sb,float  sa, \
                 float  sh, float op )
        {
            ambient[0]=ar; ambient[1]=ag; ambient[2]=ab; ambient[3]=aa;
            diffuse[0]=dr; diffuse[1]=dg; diffuse[2]=db; diffuse[3]=da;
            specular[0]=sr; specular[1]=sg; specular[2]=sb; specular[3]=sa;
            shininess=sh; opacity=op;
        }
        material(){

        }

        GLfloat ambient[4];
        GLfloat diffuse[4];
        GLfloat specular[4];
        GLfloat shininess;
        GLfloat opacity;
    };
    struct glyph{

        glyph(float  tip_r,float  tip_h,float  base_r, \
                 float  sc, int n_faces, float op )
        {
            tip[0]=tip_r; tip[1]=tip_h; base_radius=base_r;
            scale_factor=sc; N_faces= n_faces; opacity = op;
        }
        glyph(){

        }

        GLfloat tip[2]; //radius, height
        GLfloat base_radius; //radius
        GLfloat scale_factor;
        GLfloat opacity;
        GLint N_faces;


    };

    struct colour_table{
        colour_table( float r1, float r2, float r3, float r4, float r5,float r6, \
                      float g1, float g2, float g3, float g4, float g5,float g6, \
                      float b1, float b2, float b3, float b4, float b5,float b6, \
                      float a1, float a2, float a3, float a4, float a5, float a6, \
                      float sc1, float sc2, float sc3, float sc4, float sc5, float sc6, \
                      float low_clamp_r, float low_clamp_g, float low_clamp_b ,  float low_clamp_a)
        {
            r_lut[0]=r1; r_lut[1]=r2; r_lut[2]=r3; r_lut[3]=r4; r_lut[4]=r5; r_lut[5]=r6; \
            g_lut[0]=g1; g_lut[1]=g2; g_lut[2]=g3; g_lut[3]=g4; g_lut[4]=g5; g_lut[5]=g6; \
            b_lut[0]=b1; b_lut[1]=b2; b_lut[2]=b3; b_lut[3]=b4; b_lut[4]=b5; b_lut[5]=b6; \
                        a_lut[0]=a1; a_lut[1]=a2; a_lut[2]=a3; a_lut[3]=a4; a_lut[4]=a5; a_lut[5]=a6;\
            sc_lut[0]=sc1; sc_lut[1]=sc2; sc_lut[2]=sc3; sc_lut[3]=sc4; sc_lut[4]=sc5; sc_lut[5]=sc6;\
            low_clamp[0] = low_clamp_r; low_clamp[1] = low_clamp_g; low_clamp[2] = low_clamp_b;low_clamp[3] = low_clamp_a;
        }
                colour_table( float4 r, float2 r5, \
                                         float4 g, float2 g5, \
                                         float4 b, float2 b5, \
                                         float4 a, float2 a5, \
                                         float4 sc, float2 sc5, \
					 float low_clamp_r, float low_clamp_g, float low_clamp_b ,  float low_clamp_a)
        {
            r_lut[0]=r.x; r_lut[1]=r.y; r_lut[2]=r.z; r_lut[3]=r.w; r_lut[4]=r5.x; r_lut[5]=r5.y; \
            g_lut[0]=g.x; g_lut[1]=g.y; g_lut[2]=g.z; g_lut[3]=g.w; g_lut[4]=g5.x; g_lut[5]=g5.y;\
            b_lut[0]=b.x; b_lut[1]=b.y; b_lut[2]=b.z; b_lut[3]=b.w; b_lut[4]=b5.x; b_lut[5]=b5.y;\
            a_lut[0]=a.x; a_lut[1]=a.y; a_lut[2]=a.z; a_lut[3]=a.w; a_lut[4]=a5.x; a_lut[5]=a5.y;\
            sc_lut[0]=sc.x; sc_lut[1]=sc.y; sc_lut[2]=sc.z; sc_lut[3]=sc.w; sc_lut[4]=sc5.x; sc_lut[5]=sc5.y;  \
            low_clamp[0] = low_clamp_r; low_clamp[1] = low_clamp_g; low_clamp[2] = low_clamp_b;low_clamp[3] = low_clamp_a;
        }

		colour_table()
		{
                        r_lut[0]=0; r_lut[1]=0; r_lut[2]=0; r_lut[3]=0; r_lut[4]=0; r_lut[5]=0; \
            g_lut[0]=0; g_lut[1]=0; g_lut[2]=0; g_lut[3]=0; g_lut[4]=0; g_lut[5]=0;\
            b_lut[0]=0; b_lut[1]=0; b_lut[2]=0; b_lut[3]=0; b_lut[4]=0; b_lut[5]=0;\
                        a_lut[0]=0; a_lut[1]=0; a_lut[2]=0; a_lut[3]=0; a_lut[4]=0; a_lut[5]=0;\
            sc_lut[0]=0; sc_lut[1]=0; sc_lut[2]=0; sc_lut[3]=0; sc_lut[4]=0; sc_lut[5]=0; \
            low_clamp[0] = 0; low_clamp[1] = 0; low_clamp[2] = 0; low_clamp[3] = 0;
		}
		
                  GLfloat r_lut[6];
        GLfloat g_lut[6];
        GLfloat b_lut[6];
                GLfloat a_lut[6];
        GLfloat sc_lut[6];
        GLfloat low_clamp[6];
    };
	struct image_dims{
        image_dims ( float xs, float ys, float zs, float xd, float yd, float zd )
        {
            xsize=xs;
            ysize=ys;
            zsize=zs;
            xdim=xd;
            ydim=yd;
            zdim=zd;		
        }
		
        float xsize,ysize,zsize;
        float xdim,ydim,zdim;
    };
    struct bounds{
        bounds( GLfloat xmin_in, GLfloat xmax_in,  GLfloat ymin_in, GLfloat ymax_in, GLfloat zmin_in, GLfloat zmax_in)
        {
            xmin=xmin_in;
            xmax=xmax_in;
            ymin=ymin_in;
            ymax=ymax_in;
            zmin=zmin_in;
            zmax=zmax_in;
        }
        GLfloat xmin;
        GLfloat xmax;
        GLfloat ymin;
        GLfloat ymax;
        GLfloat zmin;
        GLfloat zmax;


    };

  struct vertexFloat
    {
        vertexFloat(float sc_in){
			x = 0;
			y = 0;
			z = 0;
			nx = 0;
			ny = 0;
			nz = 0;
			scalar=sc_in;
		}
        vertexFloat(){
			x = 0;
			y = 0;
			z = 0;
			nx = 0;
			ny = 0;
			nz = 0;
			scalar=0;
		}
        
        float x, y, z;        //Vertex
        float nx, ny, nz;     //Normal
        float scalar;
    };

}
#endif // BRIVIEW_STRUCTS_H
