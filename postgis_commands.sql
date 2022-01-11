
ST_Centroid (g1) /*Returns the geometric center of a geometry.*/

ST_ClosestPoint (g1, g2) /*Returns the 2-dimensional point on g1
that is closest to g2. This is the first point of the shortest line.*/

ST_Containsmm (geomA, geomB) /*Returns true if and only if no
points of B lie in the exterior of A, and at least one point of the interior
of B lies in the interior of A.*/

ST_ContainsProperly (geomA, geomB) /*Returns true if B intersects
the interior of A but not the boundary (or exterior). A does not contain
properly itself, but does contain itself.*/

ST_Crosses /*(g1, g2) Returns TRUE if the supplied geometries
have some, but not all, interior points in common.*/

ST_Disjoint (A, B) /*Returns TRUE if the Geometries do not
"spatially intersect" - if they do not share any space together.*/

ST_Distancemm (A, B) /*geometry type Returns the 2-dimensional
cartesian minimum distance (based on spatial ref) between two
geometries in projected units. For geography type defaults to return
spheroidal minimum distance between two geographies in meters.*/

ST_HausdorffDistance (A, B) /*Returns the Hausdorff distance between two
geometries. Basically a measure of how similar or dissimilar 2
geometries are. Units are in the units of the spatial reference system of
the geometries.*/

ST_MaxDistance (g1, g2) /*Returns the 2-dimensional largest
distance between two geometries in projected units.*/

ST_DFullyWithin (g1, g2, distance) /*Returns true if all of the
geometries are within the specified distance of one another*/

ST_DWithinG /*Returns true if the geometries are within the specified
distance of one another. For geometry units are in those of spatial
reference and For geography units are in meters and measurement is
defaulted to use_spheroid=true (measure around spheroid), for faster
check, use_spheroid=false to measure along sphere.*/

ST_Equalsmm (A, B) /*Returns true if the given geometries represent
the same geometry. Directionality is ignored.*/

ST_HasArc3d (geomA) /*Returns true if a geometry or geometry
collection contains a circular string.*/

-----------------------------------------------------------------------------------

GeometryType(geomA) /* Returns the type of the geometry as a string. Eg:
'LINESTRING', 'POLYGON', 'MULTIPOINT', etc.*/

ST_Boundarymm(geomA) /* Returns the closure of the combinatorial
boundary of this Geometry.*/

ST_CoordDimmm(geomA) /*Return the coordinate dimension of the
ST_Geometry value.*/

ST_Dimension2(g) /*The inherent dimension of this Geometry object, which
must be less than or equal to the coordinate dimension.*/

ST_EndPointmm(g) /*Returns the last point of a LINESTRING geometry as a
POINT.*/

ST_Envelopemm(g1) /* Returns a geometry representing the double precision
(float8) bounding box of the supplied geometry.*/

ST_ExteriorRingmm(a_polygon) /*Returns a line string representing the
exterior ring of the POLYGON geometry. Return NULL if the geometry is not a
polygon. Will not work with MULTIPOLYGON*/

ST_GeometryN(geomA, n) /*Return the 1-based Nth geometry if the
geometry is a GEOMETRYCOLLECTION, (MULTI)POINT, (MULTI)LINESTRING,
MULTICURVE or (MULTI)POLYGON, POLYHEDRALSURFACE Otherwise, return
NULL.*/

ST_GeometryType(g1) /*Return the geometry type of the ST_Geometry
value.*/

ST_NumPatches1(g1) /*Return the number of faces on a Polyhedral
Surface. Will return null for non-polyhedral geometries.*/

ST_NumPointsmm (g1) /*Return the number of points in an ST_LineString or
ST_CircularString value*/

GeometryType2(geomA) /*Returns the type of the geometry as a string. Eg:
'LINESTRING', 'POLYGON', 'MULTIPOINT', etc.*/

ST_Boundarymm(geomA) /*Returns the closure of the combinatorial
boundary of this Geometry.*/

ST_CoordDimmm(geomA) /* Return the coordinate dimension of the
ST_Geometry value.*/

ST_Dimension2(g) /*The inherent dimension of this Geometry object, which
must be less than or equal to the coordinate dimension.*/

ST_EndPointmm(g) /*Returns the last point of a LINESTRING geometry as a
POINT.*/

ST_Envelopemm(g1) /*Returns a geometry representing the double precision
(float8) bounding box of the supplied geometry.*/

ST_ExteriorRingmm(a_polygon) /*Returns a line string representing the
exterior ring of the POLYGON geometry. Return NULL if the geometry is not a
polygon. Will not work with MULTIPOLYGON*/

ST_GeometryN2(geomA, n) /*Return the 1-based Nth geometry if the
geometry is a GEOMETRYCOLLECTION, (MULTI)POINT, (MULTI)LINESTRING,
MULTICURVE or (MULTI)POLYGON, POLYHEDRALSURFACE Otherwise, return
NULL.*/

ST_GeometryType2(g1) /*Return the geometry type of the ST_Geometry
value.*/

ST_NumPatches1(g1) /*Return the number of faces on a Polyhedral
Surface. Will return null for non-polyhedral geometries.*/

ST_NumPointsmm (g1) /*Return the number of points in an ST_LineString or
ST_CircularString value.*/

--------------------------------------------
ST_Intersects /*Returns TRUE if the Geometries/Geography
"spatially intersect in 2D" - (share any portion of space) and FALSE if
they don't (they are Disjoint). For geography -- tolerance is 0.00001
meters (so any points that close are considered to intersect)*/

ST_Length /*Returns the 2d length of the geometry if it is a
linestring or multilinestring. geometry are in units of spatial reference
and geography are in meters (default spheroid)*/

ST_Length2D (a_2dlinestring) /*Returns the 2-dimensional length of
the geometry if it is a linestring or multi-linestring. This is an alias for
ST_Length*/

ST_LongestLine (g1, g2) /*Returns the 2-dimensional longest line
points of two geometries. The function will only return the first longest
line if more than one, that the function finds. The line returned will
always start in g1 and end in g2. The length of the line this function
returns will always be the same as st_maxdistance returns for g1 and
g2.*/

ST_OrderingEqualsmm (A, B) /*Returns true if the given geometries
represent the same geometry and points are in the same directional
order.*/

ST_Overlapsmm (A, B) /*Returns TRUE if the Geometries share
space, are of the same dimension, but are not completely contained by
each other.*/

ST_Perimeter /*Return the length measurement of the boundary
of an ST_Surface or ST_MultiSurface geometry or geography.
(Polygon, Multipolygon). geometry measurement is in units of spatial
reference and geography is in meters.*/

ST_Perimeter2D (geomA) /*Returns the 2-dimensional perimeter of
the geometry, if it is a polygon or multi-polygon. This is currently an
alias for ST_Perimeter.*/

ST_3DPerimeter3d (geomA) /*Returns the 3-dimensional perimeter of
the geometry, if it is a polygon or multi-polygon.*/

ST_PointOnSurfacemm(g1) /*Returns a POINT guaranteed to lie
on the surface.*/

ST_Project1(g1, distance, azimuth) /*Returns a POINT projected
from a start point using a distance in meters and bearing (azimuth) in
radians.*/

ST_Relate2 /*Returns true if this Geometry is spatially related to
anotherGeometry, by testing for intersections between the Interior,
Boundary and Exterior of the two geometries as specified by the
values in the intersectionMatrixPattern. If no intersectionMatrixPattern
is passed in, then returns the maximum intersectionMatrixPattern that
relates the 2 geometries.*/

ST_RelateMatch1(intersectionMatrix, intersectionMatrixPattern)
/*Returns true if intersectionMattrixPattern1 implies
intersectionMatrixPattern2*/

ST_ShortestLine (g1, g2) /*Returns the 2-dimensional shortest line
between two geometries*/



ST_SRID(g1) /*Returns the spatial reference identifier for the ST_Geometry
as defined in spatial_ref_sys table.
ST_StartPointmm 3d (geomA) Returns the first point of a LINESTRING
geometry as a POINT.*/

ST_Summary /*Returns a text summary of the contents of the geometry.*/

ST_Xmm(a_point) /*Return the X coordinate of the point, or NULL if not
available. Input must be a point.*/

ST_Y(a_point) /*Return the Y coordinate of the point, or NULL if not
available. Input must be a point.*/


ST_Zmflag(geomA) /*Returns ZM (dimension semantic) flag of the geometries
as a small int. Values are: 0=2d, 1=3dm, 2=3dz, 3=4d.*/

ST_ZMin3d (aGeomorBox2DorBox3D) /*Returns Z minima of a bounding box 2d
or 3d or a geometry.
1. linestring, point
2. linestring, point, position
1. geomA, a, b, c, d, e, f, g, h, i, xoff, yoff, zoff
2. geomA, a, b, d, e, xoff, yoff */


ST_ForceRHR2 /*Forces the orientation of the vertices in a polygon to
follow the Right-Hand-Rule.*/

ST_Multi (g1) /*Returns the geometry as a MULTI* geometry. If the geometry is
already a MULTI*, it is returned unchanged.*/

ST_RemovePoint(linestring, offset) /*Removes point from a linestring. Offset
is 0-based.*/

ST_Reverse (g1) /*Returns the geometry with vertex order reversed.*/

ST_Touches(g1, g2) /*Returns TRUE if the geometries have at
least one point in common, but their interiors do not intersect.*/

ST_Within(A, B) /*Returns true if the geometry A is completely
inside geometry B */

ST_Buffermm(T) /*For geometry: Returns a geometry that represents all
points whose distance from this Geometry is less than or equal to distance.
Calculations are in the Spatial Reference System of this Geometry. For
geography: Uses a planar transform wrapper. Introduced in 1.5 support for
different end cap and mitre settings to control shape. buffer_style options:
quad_segs=#,endcap=round|flat|square,join=round|mitre|bevel,mitre_limit=
1. g1, radius_of_buffer
2. g1, radius_of_buffer, num_seg_quarter_circle
3. g1, radius_of_buffer, buffer_style_parameters
4. g1, radius_of_buffer_in_meters
*/
ST_Collect  /*Return a specified ST_Geometry value from a collection of
other geometries.*/

ST_ConvexHull(geomA) /*The convex hull of a geometry represents
the minimum convex geometry that encloses all geometries within the set.
ST_CurveToLinemm 3d Converts a CIRCULARSTRING/CURVEDPOLYGON
to a LINESTRING/POLYGON*/

ST_Differencemm(geomA, geomB) /*Returns a geometry that represents
that part of geometry A that does not intersect with geometry B.*/

ST_Dump(g1) /*Returns a set of geometry_dump (geom,path) rows, that
make up a geometry g1.*/

ST_DumpPoint(geom) /*Returns a set of geometry_dump (geom,path)
rows of all points that make up a geometry.*/

ST_FlipCoordinates(geom) /*Returns a version of the given geometry
with X and Y axis flipped. Useful for people who have built latitude/longitude
features and need to fix them.*/

ST_RemoveRepeatedPoints(geom) /*Returns a version of the given
geometry with duplicated points removed*/


ST_Transform(g1, srid) /*Returns a new geometry with its coordinates
transformed to the SRID referenced by the integer parameter.*/

ST_AsBinary /*Return the Well-Known Binary (WKB) representation of
the geometry/geography without SRID meta data.*/

ST_AsEWKB2 /*Return the Well-Known Binary (WKB) representation of the
geometry with SRID meta data.*/

ST_AsEWKT2 /*Return the Well-Known Text (WKT) representation of the
geometry with SRID meta data.*/

ST_AsGeoJSONG /*Return the geometry as a GeoJSON element.
1. geom, maxdecimaldigits=15, options=0
2. geog, maxdecimaldigits=15, options=0
3. gj_version, geom, maxdecimaldigits=15, options=0
4. gj_version, geog, maxdecimaldigits=15, options=0*/

st_asgml
/*1. geom, maxdecimaldigits=15, options=0
2. geog, maxdecimaldigits=15, options=0
version, geom, maxdecimaldigits=15, options=0,
nprefix=null
3.
version, geog, maxdecimaldigits=15, options=0,
nprefix=null*/

ST_Simplify (geomA, tolerance) /*Returns a "simplified" version of the given
geometry using the Douglas-Peucker algorithm.*/

ST_SimplifyPreserveTopology (geomA, tolerance) /*Returns a "simplified"
version of the given geometry using the Douglas-Peucker algorithm. Will
avoid creating derived geometries (polygons in particular) that are invalid.
ST_Split1
 (input, blade) Returns a collection of geometries resulting by
splitting a geometry.*/

ST_SymDifferencemm(geomA, geomB) /*Returns a geometry that
represents the portions of A and B that do not intersect. It is called a
symmetric difference because ST_SymDifference(A,B) =
ST_SymDifference(B,A).*/

ST_Union  /*Returns a geometry that represents the point set union of the
Geometries.*/

ST_UnaryUnion(geom) /*Like ST_Union, but working at the
geometry component level.*/


ST_AsKML /*Return the geometry as a KML element. Several variants.
Default version=2, default precision=15*/

ST_AsSVG /*Returns a Geometry in SVG path data given a geometry or
geography object.*/

ST_AsText /*Return the Well-Known Text (WKT) representation of the
geometry/geography without SRID metadata.*/

ST_AsLatLonText
/* Return the Degrees, Minutes, Seconds representation of
the given point.*/


Box2D (geomA) /*Returns a BOX2D representing the maximum
extents of the geometry.*/









SELECT 
    st_area(w.geometry) as watershed_area_sf, 
    st_area(b.geometry) as basin_area_sf,
    st_area(w.geometry)*0.11111111111111111 as watershed_area_sy, 
    st_area(b.geometry)*0.11111111111111111 as basin_area_sy,
    st_area(w.geometry)/43560 as watershed_area_ac, 
    st_area(b.geometry)/43560 as basin_area_ac
from 
    public.watershed_sarasota as w
left join 
    public.basin_network_sarasota as b
on  
    w.objectid = b.objectid
where 
    St_intersects(b.geometry,w.geometry);

SELECT *
from public.watershed_sarasota as w
right outer join public.basin_network_sarasota as b
on  St_intersects(b.geometry,w.geometry)
where St_intersects(b.geometry,w.geometry)

