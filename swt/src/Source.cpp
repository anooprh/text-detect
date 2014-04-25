#include "mex.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
using namespace std;

#define PI 3.14159265

struct Point2d {
    int x;
    int y;
    float SWT;
};

struct Point2dFloat {
    float x;
    float y;
};

struct Ray {
    Point2d p;
    Point2d q;
    std::vector<Point2d> points;
};


void strokeWidthTransform(const float * edgeImage,
    const float * gradientX,
    const float * gradientY,
    bool dark_on_light,
    float * SWTImage,
    int h, int w,
    std::vector<Ray> & rays) {
    // First pass
    float prec = .05f;
    for( int row = 0; row < h; row++ ){
        const float* ptr = edgeImage + row*w;        
        for ( int col = 0; col < w; col++ ){
            if (*ptr > 0) {
                Ray r;

                Point2d p;
                p.x = col;
                p.y = row;
                r.p = p;
                std::vector<Point2d> points;
                points.push_back(p);

                float curX = (float)col + 0.5f;
                float curY = (float)row + 0.5f;
                int curPixX = col;
                int curPixY = row;
                float G_x = gradientX[ col + row*w ];                        
                float G_y = gradientY[ col + row*w ];
                // normalize gradient
                float mag = sqrt( (G_x * G_x) + (G_y * G_y) );
                if (dark_on_light){
                    G_x = -G_x/mag;
                    G_y = -G_y/mag;
                } else {
                    G_x = G_x/mag;
                    G_y = G_y/mag;                    
                }
                while (true) {
                    curX += G_x*prec;
                    curY += G_y*prec;
                    if ((int)(floor(curX)) != curPixX || (int)(floor(curY)) != curPixY)     {    
                        curPixX = (int)(floor(curX));
                        curPixY = (int)(floor(curY));
                        // check if pixel is outside boundary of image
                        if (curPixX < 0 || (curPixX >= w) || curPixY < 0 || (curPixY >= h)) {    
                            break;
                        }
                        Point2d pnew;
                        pnew.x = curPixX;
                        pnew.y = curPixY;
                        points.push_back(pnew);

                        if ( edgeImage[ curPixY*w+ curPixX ] > 0) {
                            r.q = pnew;
                            // dot product
                            float G_xt = gradientX[ curPixY*w + curPixX ];
                            float G_yt = gradientY[ curPixY*w + curPixX ];
                            mag = sqrt( (G_xt * G_xt) + (G_yt * G_yt) );
                            if (dark_on_light){
                                G_xt = -G_xt/mag;
                                    G_yt = -G_yt/mag;
                            } else {
                                G_xt = G_xt/mag;
                                G_yt = G_yt/mag;                                
                            }

                            if (acos(G_x * -G_xt + G_y * -G_yt) < PI/2.0 ) {
                                float length = sqrt( ((float)r.q.x - (float)r.p.x)*((float)r.q.x - (float)r.p.x) + ((float)r.q.y - (float)r.p.y)*((float)r.q.y - (float)r.p.y));
                                for (std::vector<Point2d>::iterator pit = points.begin(); pit != points.end(); pit++) {
                                    float* pSWT = SWTImage +  w * pit->y + pit->x;
                                    if (*pSWT < 0) {
                                        *pSWT = length;
                                    } else {
                                        *pSWT = std::min(length, *pSWT);
                                    }
                                }
                                r.points = points;
                                rays.push_back(r);
                            }
                            break;
                        }
                    }
                }
            }
            ptr++;
        }
    }    
}


bool Point2dSort(const Point2d &lhs, const Point2d &rhs) {
    return lhs.SWT < rhs.SWT;
}

void SWTMedianFilter(float * SWTImage, int h, int w,
        std::vector<Ray> & rays, float maxWidth = -1 ) {
    for (std::vector<Ray>::iterator rit = rays.begin(); rit != rays.end(); rit++) {
        for (std::vector<Point2d>::iterator pit = rit->points.begin(); pit != rit->points.end(); pit++) {
            pit->SWT = SWTImage[ w*pit->y + pit->x ];
        }
        std::sort(rit->points.begin(), rit->points.end(), &Point2dSort);
        //std::nth_element( rit->points.begin(), rit->points.end(), rit->points.size()/2, &Point2dSort );
        float median = (rit->points[rit->points.size()/2]).SWT;
        if ( maxWidth > 0 && median >= maxWidth ) {
            median = -1;
        }
        for (std::vector<Point2d>::iterator pit = rit->points.begin(); pit != rit->points.end(); pit++) {
            SWTImage[ w*pit->y + pit->x ] = std::min(pit->SWT, median);
        }
    }    
}

typedef std::vector< std::set<int> > graph_t; // graph as a list of neighbors per node

void connComp( const graph_t& g, std::vector<int>& c, int i, int l ) {
    // starting from node i labe this conn-comp with label l
    if ( i < 0 || i > g.size() ) {
        return;
    }
    std::vector< int > stack;
    // push i
    stack.push_back(i);
    c[i] = l;
    while ( ! stack.empty() ) {
        // pop
        i = stack.back();
        stack.pop_back();
        // go over all nieghbors
        for ( std::set<int>::const_iterator it = g[i].begin(); it != g[i].end(); it++ ) {
            if ( c[*it] < 0 ) {
                stack.push_back( *it );
                c[ *it ] = l;
            }
        }
    }
}
int findNextToLabel( const graph_t& g, const vector<int>& c ) {
    for ( int i = 0 ; i < c.size(); i++ ) {
        if ( c[i] < 0 ) {
            return i;
        }
    }
    return c.size();
}

int connected_components(const graph_t& g, vector<int>& c) {
    // check for empty graph!
    if ( g.empty() ) {
        return 0;
    }
    int i = 0;
    int num_conn = 0;
    do {
        connComp( g, c, i, num_conn );
        num_conn++;
        i = findNextToLabel( g, c );
    } while ( i < g.size() );
    return num_conn;
}

std::vector< std::vector<Point2d> >
        findLegallyConnectedComponents(const float* SWTImage, int h, int w,
        std::vector<Ray> & rays) {
    std::map<int, int> Map;
    std::map<int, Point2d> revmap;
    std::vector<std::vector<Point2d> > components; // empty
    int num_vertices = 0, idx = 0;
    graph_t g;
    // Number vertices for graph.  Associate each point with number
    for( int row = 0; row < h; row++ ){        
        for (int col = 0; col < w; col++ ){
            idx = col + w * row;
            if (SWTImage[idx] > 0) {
                Map[idx] = num_vertices;
                Point2d p;
                p.x = col;
                p.y = row;
                revmap[num_vertices] = p;
                num_vertices++;
                std::set<int> empty;
                g.push_back(empty);
            }
        }
    }   
    if ( g.empty() ) {
        return components; // nothing to do with an empty graph...
    }
    for( int row = 0; row < h; row++ ){        
        for (int col = 0; col < w; col++ ){
            idx = col + w * row;
            if ( SWTImage[idx] > 0) {
                // check pixel to the right, right-down, down, left-down
                int this_pixel = Map[idx];
                float thisVal = SWTImage[idx];
                if (col+1 < w) {
                    float right = SWTImage[ w*row + col + 1 ];
                    if (right > 0 && (thisVal/right <= 3.0 || right/thisVal <= 3.0)) {
                        g[this_pixel].insert( Map[ w*row + col + 1 ] );                    
                        g[ Map[ w*row + col + 1 ] ].insert( this_pixel );
                        //boost::add_edge(this_pixel, map.at(row * SWTImage->width + col + 1), g);
                    }
                }
                if (row+1 < h) {
                    if (col+1 < w) {
                        float right_down = SWTImage[ w*(row+1) + col + 1 ];
                        if (right_down > 0 && (thisVal/right_down <= 3.0 || right_down/thisVal <= 3.0)) {
                            g[ this_pixel ].insert( Map[ w*(row+1) + col + 1 ] );
                            g[ Map[ w*(row+1) + col + 1 ] ].insert(this_pixel);                            
                            // boost::add_edge(this_pixel, map.at((row+1) * SWTImage->width + col + 1), g);
                        }
                    }
                    float down = SWTImage[ w*(row+1) + col ];
                    if (down > 0 && (thisVal/down <= 3.0 || down/thisVal <= 3.0)) {
                        g[ this_pixel ].insert( Map[ w*(row+1) + col ] );
                        g[ Map[ w*(row+1) + col ] ].insert( this_pixel );
                        //boost::add_edge(this_pixel, map.at((row+1) * SWTImage->width + col), g);
                    }
                    if (col-1 >= 0) {
                        float left_down = SWTImage[ w*(row+1) + col - 1 ];
                        if (left_down > 0 && (thisVal/left_down <= 3.0 || left_down/thisVal <= 3.0)) {
                            g[ this_pixel ].insert( Map[ w*(row+1) + col - 1 ] );
                            g[ Map[ w*(row+1) + col - 1 ] ].insert( this_pixel );
                            //boost::add_edge(this_pixel, map.at((row+1) * SWTImage->width + col - 1), g);
                        }
                    }
                }    
            }            
        }
    }

    std::vector<int> c(num_vertices, -1);    
    int num_comp = connected_components(g, c);    

    components.reserve(num_comp);
    //std::cout << "Before filtering, " << num_comp << " components and " <<     num_vertices << " vertices" << std::endl;
    for (int j = 0; j < num_comp; j++) {
        std::vector<Point2d> tmp;
        components.push_back( tmp );
    }
    for (int j = 0; j < num_vertices; j++) {
        Point2d p = revmap[j];
        (components[c[j]]).push_back(p);
    }

    return components;
}

enum {
    EIN = 0,
    GXIN,
    GYIN,
    DOLFIN,
    MAXWIN,
    NIN };

void mexFunction( int nout, mxArray* pout[], int nin, const mxArray* pin[] ) {
    //    
    // make sure images are input in transposed so that they are arranged row-major in memory
    //
    mxAssert( nin == NIN, "wrong number of inputs" );
    mxAssert( nout > 1, "only one output" );

    int h = mxGetN( pin[EIN] ); // inputs are transposed!
    int w = mxGetM( pin[EIN] );

    mxAssert( mxIsClass( pin[EIN], mxSINGLE_CLASS ) && h == mxGetN( pin[EIN] ) && w == mxGetM( pin[EIN] ), "edge map incorrect");
    mxAssert( mxIsClass( pin[GXIN], mxSINGLE_CLASS ) && h == mxGetN( pin[GXIN] ) && w == mxGetM( pin[GXIN] ), "edge map incorrect");
    mxAssert( mxIsClass( pin[GYIN], mxSINGLE_CLASS ) && h == mxGetN( pin[GYIN] ) && w == mxGetM( pin[GYIN] ), "edge map incorrect");

    const float * edgeImage = (float*) mxGetData( pin[EIN] );
    const float * gradientX = (float*) mxGetData( pin[GXIN] );
    const float * gradientY = (float*) mxGetData( pin[GYIN] );

    bool dark_on_light = mxGetScalar( pin[DOLFIN] ) != 0 ;
    float maxWidth = mxGetScalar( pin[MAXWIN] );

    // allocate output
    pout[0] = mxCreateNumericMatrix( w, h, mxSINGLE_CLASS, mxREAL );
    float * SWTImage = (float*) mxGetData( pout[0] );
    // set SWT to -1
    for ( int i = 0 ; i < w*h; i++ ) {
        SWTImage[i] = -1;
    }

    std::vector<Ray> rays;
    strokeWidthTransform ( edgeImage, gradientX, gradientY, dark_on_light, SWTImage, h, w, rays );
    SWTMedianFilter ( SWTImage, h, w, rays, maxWidth );

    // connected components
    if ( nout > 1 ) {
        // Calculate legally connect components from SWT and gradient image.
        // return type is a vector of vectors, where each outer vector is a component and
        // the inner vector contains the (y,x) of each pixel in that component.
        std::vector<std::vector<Point2d> > components = findLegallyConnectedComponents(SWTImage, h, w, rays);
        pout[1] = mxCreateNumericMatrix( w, h, mxSINGLE_CLASS, mxREAL );
        float* pComp = (float*) mxGetData( pout[1] );
        for ( int i = 0 ; i < w*h; i++ ) {
            pComp[i] = 0;
        }
        for ( int ci = 0 ; ci < components.size(); ci++ ) {
            for ( std::vector<Point2d>::iterator it = components[ci].begin() ; it != components[ci].end(); it++ ) {
                pComp[ w * it->y + it->x ] = ci + 1;
            }
        }
    }
}