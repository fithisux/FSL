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
#include <fslsurfacegl.h>
#include	<iostream>

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

using namespace std;

namespace fslsurface_name{

    // window identifier
    static int glut_win;
    //static int win=0;
    const fslSurface<float,unsigned int>* surf2render;
    static GLuint glut_vbos[2];
    static  GLint InVertex_loc, InNormal_loc, InScalar_loc;
    
    static char* textFileRead(const char *fileName) {
        char* text;
        
        if (fileName != NULL) {
            FILE *file = fopen(fileName, "rt");
            
            if (file != NULL) {
                fseek(file, 0, SEEK_END);
                int count = ftell(file);
                rewind(file);
                
                if (count > 0) {
                    text = (char*)malloc(sizeof(char) * (count + 1));
                    count = fread(text, sizeof(char), count, file);
                    text[count] = '\0';
                }
                fclose(file);
            }
        }
        return text;
    }

    
    
    
    void glutRender( const fslSurface<float, unsigned int> & surf , string name )
    {
        int argc=0;
        char* argv="";
        glutInitAndCreateWindow(&argc, &(argv), name);
        glutDisplayFunc(glutDisp);
        surf2render=&surf;
        
        glGenBuffersARB(2,glut_vbos);
        glBufferData_Vertices(*surf2render, glut_vbos[0]);
        glBufferData_Faces(*surf2render, glut_vbos[1]);
       // InVertex_loc=0;
      //  InNormal_loc=1;
     //   InScalar_loc=2;
        
        GLuint v_light_dir_map_scalars, p_light_dir_map_scalars,f_light_dir_map_scalars;
        v_light_dir_map_scalars = glCreateShader(GL_VERTEX_SHADER);
        f_light_dir_map_scalars = glCreateShader(GL_FRAGMENT_SHADER);
        cout<<"read shaders"<<endl;
         char *vs,*fs;
        vs = textFileRead(("/Users/brianpatenaude/fslsrc/briview/briview.app/Contents/MacOS/glsl_shaders/surface_dir_light_map_scalars_ctrthresh.vert"));        
        fs = textFileRead(("/Users/brianpatenaude/fslsrc/briview/briview.app/Contents/MacOS/glsl_shaders/surface_dir_light_map_scalars.frag"));
        cout<<"done reading shaders "<<endl;
        
        const char * vv_light_dir = vs;
        const char * ff_light_dir = fs;
        
        glShaderSource(v_light_dir_map_scalars, 1, &vv_light_dir,NULL);
        glShaderSource(f_light_dir_map_scalars, 1, &ff_light_dir,NULL);
        
        glCompileShader(v_light_dir_map_scalars);
        glCompileShader(f_light_dir_map_scalars);
        
        p_light_dir_map_scalars = glCreateProgram();
        
        glAttachShader(p_light_dir_map_scalars,v_light_dir_map_scalars);
        glAttachShader(p_light_dir_map_scalars,f_light_dir_map_scalars);
        
        glBindAttribLocation(p_light_dir_map_scalars,0,"InVertex");
        glBindAttribLocation(p_light_dir_map_scalars,1,"InNormal");
        glBindAttribLocation(p_light_dir_map_scalars,2,"InScalar");
        
        
        InVertex_loc = 0;//glGetAttribLocation(p_light_dir_map_scalars, "InVertex");
        InNormal_loc = 1;//glGetAttribLocation(p_light_dir_map_scalars, "InNormal");
        InScalar_loc = 2;//glGetAttribLocation(p_light_dir_map_scalars, "InScalar");
        glLinkProgram(p_light_dir_map_scalars);
        glUseProgram(p_light_dir_map_scalars);
        
        GLint loc_r_lut,loc_g_lut,loc_b_lut,loc_a_lut, loc_sc_lut;
        GLint loc_r_lut_last,loc_g_lut_last,loc_b_lut_last,loc_a_lut_last, loc_sc_lut_last;
        GLint loc_low_clamp;
        GLfloat r_lut[4],g_lut[4],b_lut[4],a_lut[4],sc_lut[4];
        GLfloat r_lut_last[2],g_lut_last[2],b_lut_last[2],a_lut_last[2], sc_lut_last[2];

        
        loc_r_lut = glGetUniformLocation(p_light_dir_map_scalars,"r_lut");
        loc_g_lut = glGetUniformLocation(p_light_dir_map_scalars,"g_lut");
        loc_b_lut = glGetUniformLocation(p_light_dir_map_scalars,"b_lut");
        loc_a_lut = glGetUniformLocation(p_light_dir_map_scalars,"a_lut");
        
        loc_sc_lut = glGetUniformLocation(p_light_dir_map_scalars,"sc_lut");
        loc_r_lut_last = glGetUniformLocation(p_light_dir_map_scalars,"r_lut_last");
        loc_g_lut_last = glGetUniformLocation(p_light_dir_map_scalars,"g_lut_last");
        loc_b_lut_last = glGetUniformLocation(p_light_dir_map_scalars,"b_lut_last");
        loc_a_lut_last = glGetUniformLocation(p_light_dir_map_scalars,"a_lut_last");
        
        loc_sc_lut_last = glGetUniformLocation(p_light_dir_map_scalars,"sc_lut_last");
        loc_low_clamp = glGetUniformLocation(p_light_dir_map_scalars,"lower_clamp");
     
        r_lut[0]=1.0;
        r_lut[1]=1.0;
        r_lut[2]=0.0;
        r_lut[3]=0.0;
        r_lut_last[0]=0.0;
        r_lut_last[1]=0.0;
        
        g_lut[0]=0.0;
        g_lut[1]=1.0;
        g_lut[2]=1.0;
        g_lut[3]=1.0;
        g_lut_last[0]=0.0;
        g_lut_last[1]=0.0;
        
        b_lut[0]=0.0;
        b_lut[1]=0.0;
        b_lut[2]=0.0;
        b_lut[3]=1.0;
        b_lut_last[0]=1.0;
        b_lut_last[1]=1.0;
        
        a_lut[0]=1.0;
        a_lut[1]=1.0;
        a_lut[2]=1.0;
        a_lut[3]=1.0;
        a_lut_last[0]=1.0;
        a_lut_last[1]=1.0;
        
        sc_lut[0]=0.0;
        sc_lut[1]=1.0/5;
        sc_lut[2]=2.0/5;
        sc_lut[3]=3.0/5;
        sc_lut_last[0]=4.0/5;
        sc_lut_last[1]=5.0/5;
        
        glUniform4fv(loc_r_lut,1,r_lut);
        glUniform4fv(loc_g_lut,1,g_lut);
        glUniform4fv(loc_b_lut,1,b_lut);
        glUniform4fv(loc_a_lut,1,a_lut);
        glUniform4fv(loc_sc_lut,1,sc_lut);
        
        glUniform2fv(loc_r_lut_last,1,r_lut_last);
        glUniform2fv(loc_g_lut_last,1,g_lut_last);
        glUniform2fv(loc_b_lut_last,1,b_lut_last);
        glUniform2fv(loc_a_lut_last,1,a_lut_last);
        glUniform2fv(loc_sc_lut_last,1,sc_lut_last);
        
        glUniform4f(loc_low_clamp,r_lut[0],g_lut[0],b_lut[0],a_lut[0]);
      //cout<<"locations "<<InVertex_loc<<" "<<InNormal_loc<<" "<<InScalar_loc<<endl;
     //   surfaces->setVertexAttribLocs(InVertex_loc,InNormal_loc, InScalar_loc);
        
        
        cout<<surf2render->getNumberOfFaces()<<endl;
      //  glutIdleFunc(glutIdle);
        // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glutKeyboardFunc(keyb);
        
        //  glBegin(GL_TRIANGLES);
        //glVertex3f(-0.5,-0.5,0.0);
        //glVertex3f(0.5,0.0,0.0);
        //glVertex3f(0.0,0.5,0.0);
        //glEnd();
        glutMainLoop();
     //   glutDestroyWindow(win);
    }
    void glutInitAndCreateWindow(int *argc, char **argv, const std::string & name ){
        glutInit(argc, argv);
        glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
        glutInitWindowPosition(100,100);
        glutInitWindowSize(320,320);
        glut_win = glutCreateWindow(name.c_str());
    }
    void glutDisp( void )
    {
        cout<<"glutDisp "<<endl;
        glClearColor(0.0,0.0,0.0,1.0);

        glClear(GL_COLOR_BUFFER_BIT);
        render<float,unsigned int>(*surf2render, InVertex_loc, InNormal_loc, InScalar_loc,glut_vbos[0] ,glut_vbos[1] );

        glutSwapBuffers();

        
    }
  //  void glutIdle( void )
   // {
       // cout<<"glutIdle "<<endl;
     //   glClearColor(1.0,0.0,0.0,1.0);
        
       // glClear(GL_COLOR_BUFFER_BIT);
        //glutSwapBuffers();
        
   // }
    
    void keyb(unsigned char key, int x, int y){
        cout << "Pressed key " << key << " on coordinates (" << x << "," << y << ")";
        cout << endl;
        if(key == 'e'){
            cout << "Got e,so quitting " << endl;
            glutDestroyWindow(glut_win);
            glDeleteBuffersARB(2,glut_vbos);

            exit(0);
        }
    }
    
    template<class T, class T2>
	void render( const fslSurface<T,T2> & surf, const GLint & vertexLoc , const GLint &  normalLoc, const GLint & scalarLoc, const GLuint & vbo_verts, const GLuint & vbo_ele_array ) 
	{
        //	//cout<<"renderin "<<" "<<vertexLoc<<" "<<normalLoc<<" "<<scalarLoc<<endl;
        //cout<<"Render "<<surf->getNumberOfVertices()<<" "<<surf->getNumberOfFaces()<<endl;
        //	cout<<"render "<<endl;
		glBindBuffer(GL_ARRAY_BUFFER, vbo_verts );         // for vertex coordinates
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_ele_array );
		
		//cout<<"render2 "<<endl;
#ifdef __linux
		glEnableVertexAttribArrayARB(vertexLoc);
		glEnableVertexAttribArrayARB(normalLoc);
		glEnableVertexAttribArrayARB(scalarLoc);
#else
		glEnableVertexAttribArray(vertexLoc);
		glEnableVertexAttribArray(normalLoc);
		glEnableVertexAttribArray(scalarLoc);
#endif
		
		//cout<<"render3 "<<endl;
#ifdef __linux
		
		glVertexAttribPointerARB(vertexLoc, 3, GL_FLOAT, GL_FALSE, sizeof(vertex<T>), 0);
		glVertexAttribPointerARB(normalLoc, 3, GL_FLOAT, GL_FALSE, sizeof(vertex<T>), BUFFER_OFFSET(3*sizeof(float)));
		glVertexAttribPointerARB(scalarLoc, 1, GL_FLOAT, GL_FALSE, sizeof(vertex<T>), BUFFER_OFFSET(6*sizeof(float)));
#else
		glVertexAttribPointer(vertexLoc, 3, GL_FLOAT, GL_FALSE, sizeof(vertex<T>), 0);
		glVertexAttribPointer(normalLoc, 3, GL_FLOAT, GL_FALSE, sizeof(vertex<T>), BUFFER_OFFSET(3*sizeof(float)));
		glVertexAttribPointer(scalarLoc, 1, GL_FLOAT, GL_FALSE, sizeof(vertex<T>), BUFFER_OFFSET(6*sizeof(float)));
		
#endif
		
		//cout<<"render4 "<<surf->N_triangles<<endl;
		//	glDrawElements(GL_TRIANGLES,1*3,GL_UNSIGNED_INT,0);
		
		glDrawElements(GL_TRIANGLES,surf.N_triangles*3,GL_UNSIGNED_INT,0);
		//cout<<"render5 "<<surf->N_triangles<<endl;
		
#ifdef __linux
		glDisableVertexAttribArrayARB(vertexLoc);
		glDisableVertexAttribArrayARB(normalLoc);
		glDisableVertexAttribArrayARB(scalarLoc);
		
#else
		glDisableVertexAttribArray(vertexLoc);
		glDisableVertexAttribArray(normalLoc);
		glDisableVertexAttribArray(scalarLoc);
		
#endif
		
		glBindBuffer(GL_ARRAY_BUFFER, 0 );         // for vertex coordinates
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0 );
	}
	

	template<class T, class T2>
	void glBufferData_Vertices( const fslSurface<T,T2> & surf, const GLuint & vbo )  
	{
		////cout<<"buffer data "<<vertices[0].nx<<" "<<vertices[0].ny<<" "<<vertices[0].nz<<endl;
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, surf.N_vertices * sizeof(vertex<T>), &(surf.vertices[0]), GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		
	}

	
	template<class T, class T2>
	void glBufferSubData_Vertices( const fslSurface<T,T2> & surf, const GLuint & vbo )  
	{
		//cout<<"buffer data "<<surf->N_vertices<<" "<<surf->vertices.size()<<endl;;//vertices[0].nx<<" "<<vertices[0].ny<<" "<<vertices[0].nz<<endl;
		
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0,  sizeof(vertex<T>)*surf.N_vertices, &surf.vertices[0]);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		
	}	
    
    template<class T,class T2>
    void depthSortTriangles( const fslSurface<T,T2> & surf, const GLuint & vbo_vert, const GLuint & vbo_tris ){

        vertex<T>* verts = new vertex<T>[surf.N_vertices];
        T2* tris = new T2[surf.N_triangles*3];

        glBindBuffer(GL_ARRAY_BUFFER, vbo_vert);
        glGetBufferSubData(	GL_ARRAY_BUFFER,0,surf.N_vertices * sizeof(vertex<T>),verts);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
   /*     
        for (unsigned int i_v = 0; i_v != surf->N_vertices; ++i_v)
        {
            cout<<"vert "<<verts[i_v].x<<" "<<verts[i_v].y<<" "<<verts[i_v].z<<endl;
        }
        
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_tris);
        glGetBufferSubData(	GL_ELEMENT_ARRAY_BUFFER,0,surf->N_triangles *3* sizeof(T2),tris);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        
        for (unsigned int i_v = 0; i_v != surf->N_triangles*3; i_v+=3)
        {
            cout<<"tris "<<tris[i_v]<<" "<<tris[i_v+1]<<" "<<tris[i_v+2]<<endl;
        }
*/
        cout<<"depth sort temp"<<endl;
        
        delete verts;

    }

	
	
	template<class T, class T2>
	void glBufferData_Faces( const fslSurface<T,T2> & surf, const GLuint & vbo )  
	{		
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,vbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, surf.N_triangles * 3 * sizeof(T2), &(surf.faces[0]), GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		
	}

	template<class T, class T2>
	void glBufferSubData_Faces( const fslSurface<T,T2> & surf, const GLuint & vbo )  
	{		
	  ////cout<<"buffer sub data "<<surf->N_triangles * 3<<" "<<surf->faces.size()<<endl;
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo);
		//	unsigned int f=1;
		//glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,8,&f);
		glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, surf.N_triangles * 3 * sizeof(T2),&surf.faces[0]);
		
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	
	
	template void glBufferData_Vertices<float, unsigned int>( const fslSurface<float, unsigned int> & surf, const GLuint & vbo );  
	template void glBufferSubData_Vertices<float, unsigned int>( const fslSurface<float, unsigned int> & surf, const GLuint & vbo )  ;
	template 	void glBufferData_Faces<float, unsigned int>( const fslSurface<float, unsigned int> & surf, const GLuint & vbo )  ;
	template void glBufferSubData_Faces<float, unsigned int>( const fslSurface<float, unsigned int> & surf, const GLuint & vbo )  ;
	template void render<float, unsigned int>( const fslSurface<float, unsigned int> & surf, const GLint & vertexLoc , const GLint &  normalLoc, const GLint & scalarLoc, const GLuint & vbo_verts, const GLuint & vbo_ele_array ) ;
    template void depthSortTriangles<float, unsigned int>(const fslSurface<float, unsigned int> & surf, const GLuint & vbo_verts,const GLuint & vbo_tris);

	
	

}
