package amenon;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.atan;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;

import javax.xml.parsers.ParserConfigurationException;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateArrays;
import org.locationtech.jts.geom.CoordinateSequence;
import org.locationtech.jts.geom.CoordinateSequenceFilter;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.IntersectionMatrix;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.MultiPolygon;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKBReader;
import org.locationtech.jts.io.WKBWriter;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.io.WKTWriter;
import org.locationtech.jts.io.geojson.GeoJsonWriter;
import org.locationtech.jts.io.gml2.GMLReader;
import org.locationtech.jts.io.gml2.GMLWriter;
import org.locationtech.jts.io.kml.KMLWriter;

import org.locationtech.jts.io.geojson.GeoJsonReader;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.locationtech.jts.operation.overlay.snap.GeometrySnapper;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.operation.valid.IsValidOp;
import org.locationtech.jts.operation.valid.TopologyValidationError;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.locationtech.jts.util.GeometricShapeFactory;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.algorithm.MinimumBoundingCircle;
import org.locationtech.jts.algorithm.MinimumDiameter;
import org.locationtech.jts.densify.Densifier;
import org.locationtech.proj4j.CRSFactory;
import org.locationtech.proj4j.CoordinateReferenceSystem;
import org.locationtech.proj4j.CoordinateTransform;
import org.locationtech.proj4j.CoordinateTransformFactory;
import org.locationtech.proj4j.ProjCoordinate;
import org.xml.sax.SAXException;

public class GeoFunctions {
	
	private static final int NO_SRID = 0;
	private static final SpatialReference SPATIAL_REFERENCE = SpatialReference.create(4326);

	GeoFunctions() {
	}
	
	/*
	 * public static void main(String args[]) {
	 * 
	 * Coordinate[] coordinates = new Coordinate[]{new Coordinate(0, 0), new
	 * Coordinate(10, 10), new Coordinate(20, 20)}; // use the default factory,
	 * which gives full double-precision Geometry g2 = new
	 * GeometryFactory().createLineString(coordinates);
	 * System.out.println("Geometry 2: " + g2);
	 * 
	 * GeometryEngine geng= new GeometryEngine();
	 * 
	 * GeoFunctions gf= new GeoFunctions(); // String s=geng.googlemaplink(g2, "m",
	 * 19);
	 * 
	 * // System.out.println("geo to google is--->"+s);
	 * 
	 * String wkt =
	 * "MULTIPOLYGON (((355815.5 6689436.9, 355816.1 6689439.3, 355817.3 6689441.4, 355818.9 6689443.3, 355820.7 6689444.5, 355822.9 6689445.3, 355827.7 6689445.7, 355831.9 6689458.7, 355829.8 6689459.9, 355828 6689461.2, 355826.6 6689463.5, 355826.2 6689466.1, 355826.7 6689468.2, 355827.5 6689469.6, 355829.1 6689471.1, 355831.5 6689472.1, 355834.3 6689472.1, 355836.5 6689471.4, 355838.9 6689477.4, 355838 6689477.6, 355835.9 6689478.6, 355834.3 6689480.2, 355833.8 6689481.2, 355833.2 6689482.7, 355833.1 6689484.8, 355833.7 6689486.9, 355835 6689488.6, 355837 6689490, 355839.7 6689490.5, 355841.8 6689490.2, 355843.2 6689489.4, 355845.4 6689488.6, 355848.2 6689495.7, 355852.5 6689494.1, 355850 6689488.2, 355849.5 6689488.3, 355849 6689487, 355847.1 6689482.3, 355854.6 6689479.4, 355855.3 6689480.9, 355856.3 6689480.5, 355842.3 6689445.5, 355844.6 6689443.6, 355845.2 6689441.2, 355846.1 6689440.9, 355847.2 6689441.2, 355847 6689441.8, 355847.6 6689442, 355847.8 6689441.2, 355849.5 6689440.5, 355850.2 6689440.7, 355850.5 6689440.1, 355850.1 6689440.1, 355850.7 6689437.6, 355851.1 6689437.6, 355851.1 6689437.1, 355850.4 6689437, 355859 6689407.7, 355862.3 6689396.8, 355866 6689397.7, 355866.4 6689396.9, 355907 6689416.9, 355904.2 6689423.3, 355930.6 6689436.3, 355935 6689430.8, 355940.2 6689434.5, 355940.4 6689435.3, 355939.2 6689435.8, 355955.1 6689485.7, 355958.1 6689485.1, 355958.1 6689486.5, 355958.9 6689487.5, 355960.5 6689495.3, 355957.3 6689494.8, 355956.8 6689493.3, 355931.7 6689498.4, 355926.9 6689498.8, 355903.4 6689502.1, 355904.3 6689499.5, 355901.4 6689494.7, 355901.7 6689490.7, 355897.4 6689483.6, 355889.4 6689488.1, 355889.7 6689488.6, 355888 6689490.2, 355887.5 6689490, 355887.2 6689490.6, 355887.6 6689490.7, 355888.1 6689493.3, 355887.7 6689493.5, 355887.8 6689494, 355887 6689494.3, 355887.3 6689495.3, 355889.5 6689494.8, 355890.5 6689495.3, 355891.9 6689500.4, 355889.1 6689500.9, 355889.4 6689502.6, 355889.9 6689502.5, 355891.8 6689510.4, 355889.8 6689510.6, 355889.9 6689514.8, 355886.9 6689515.7, 355885.6 6689522.9, 355907 6689528.7, 355903.8 6689511.1, 355904.9 6689511, 355906 6689510.5, 355906.1 6689509.3, 355905.9 6689508.5, 355933 6689506.3, 355954.3 6689500.6, 355957 6689502.7, 355975.4 6689503.8, 355978.9 6689503.7, 355980.8 6689503.4, 355983.6 6689502.4, 355986.1 6689500.7, 355987.5 6689499.2, 355988.6 6689497.5, 355989.7 6689494.8, 355990.2 6689491.8, 355990 6689489.8, 355989.3 6689487.2, 355988.4 6689485.4, 355987.3 6689483.9, 355984.8 6689481.7, 355982.2 6689480.4, 355977.1 6689480.1, 355958.7 6689481.6, 355944.4 6689437.9, 355947.3 6689431.9, 355947.5 6689430.7, 355947.5 6689428.2, 355946.5 6689425.8, 355945.9 6689425, 355943.8 6689423.4, 355942.3 6689422.8, 355940.5 6689422.5, 355938.6 6689422.7, 355937 6689423.4, 355935.5 6689424.4, 355934.4 6689425.8, 355933.6 6689427.4, 355910.2 6689415.8, 355908.1 6689414.7, 355860.1 6689390.1, 355859 6689386, 355858.4 6689384.6, 355857.2 6689382.9, 355856 6689382, 355854.7 6689381.3, 355851.7 6689380.7, 355848.7 6689381.5, 355847.9 6689382, 355846.2 6689383.6, 355845.5 6689385, 355845.1 6689387.1, 355845.3 6689388.9, 355845.7 6689390.4, 355846.4 6689391.9, 355847.6 6689393.2, 355849 6689394.1, 355850 6689394.5, 355846.6 6689401.7, 355846.2 6689402.6, 355838 6689429.1, 355831.5 6689427.5, 355826.7 6689426.4, 355824.2 6689426.4, 355821.7 6689426.9, 355819.5 6689428.1, 355817.8 6689429.7, 355816.5 6689431.6, 355815.7 6689433.9, 355815.5 6689436.9)), ((355932.4 6689440.8, 355915.2 6689446.6, 355929 6689490.4, 355946.4 6689484.8, 355932.4 6689440.8)))"
	 * ; String wkt2="POINT (1620181.9310988 -2673241.00461511)"; Geom g =
	 * GeoFunctions.ST_GeomFromText(wkt, 2154); Geom g1 =
	 * GeoFunctions.ST_GeomFromText(wkt2, 3577);
	 * 
	 * String result = GeoFunctions.ST_GoogleMapLink(g1, "m", 19);
	 * 
	 * System.out.println("geo to google from wkt is--->"+result);
	 * 
	 * }
	 */
	
	
	public static Geom ST_GeomFromText(String s) {
		return ST_GeomFromText(s, NO_SRID);
	}
	
	public static Geom ST_GeomFromText(String s, int srid) {
		final Geometry g = GeometryEngine.geometryFromWkt(s, GeometryType.GEOMETRY);
		return bind(g, srid);
	}
	
	public static String ST_GoogleMapLink(Geom geom, String layerType, int zoom) {
		String result = GeometryEngine.googlemaplink(geom.g(), layerType, zoom);
		return result;
	}
	
	public static List<Double> ST_Transform(Geom geom, int srid) {
		List<Double> latlong= new ArrayList<Double>();
		Geometry convertedGeometry = GeometryEngine.transform(geom.g(), srid);
		 Coordinate c1 = convertedGeometry.getEnvelopeInternal().centre();
		 Double latit=c1.y;
		 Double longitude=c1.x;
		 latlong.add(latit);
		 latlong.add(longitude);
		return latlong;
	}
	
	static class GeometryEngine {
		
		private static final WKTReader wkr = new WKTReader();
		private static UnsupportedOperationException err() {
			return new UnsupportedOperationException();
		}
		
		
		
		//converts wkt to geometry
		
		public static Geometry geometryFromWkt(String wkt, GeometryType geometrytype) {
			if (wkt == null) {
				return null;
			}
			Geometry geom = null;
			try {
				geom = wkr.read(wkt);
			} catch (ParseException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			switch (geometrytype) {
			case GEOMETRY:

				break;
			case POINT:
				if (!geom.getGeometryType().equalsIgnoreCase("point")) {
					throw err();
				}
				break;
			case LINESTRING:
				if (!geom.getGeometryType().equalsIgnoreCase("linestring")) {
					throw err();
				}
				break;
			case POLYGON:
				if (!geom.getGeometryType().equalsIgnoreCase("polygon")) {
					throw err();
				}
				break;
			case MULTIPOINT:
				if (!geom.getGeometryType().equalsIgnoreCase("MULTIPOINT")) {
					throw err();
				}
				break;

			case MULTILINESTRING:
				if (!geom.getGeometryType().equalsIgnoreCase("MULTILINESTRING")) {
					throw err();
				}
				break;

			case MULTIPOLYGON:
				if (!geom.getGeometryType().equalsIgnoreCase("MULTIPOLYGON")) {
					throw err();
				}
				break;

			default:
				throw err();

			}
			return geom;

		}
		
		
		public static String googlemaplink(Geometry geom, String layer, int zoomlevel) {
			String result = null;

			if (geom == null) {
				return null;
			}
			if (layer == null) {
				layer = "m";
			} else if (!layer.equals("m") && !layer.equals("k") && !layer.equals("h") && !layer.equals("p")) {
				return null;
			}
			if (String.valueOf(zoomlevel) == null) {
				zoomlevel = 19;
			} else if (zoomlevel > 19 || zoomlevel < 1) {
				return null;
			}
			Geometry g = null;
			if (geom.getSRID() != 4326) {
				g = transform(geom, 4326);
			} else {
				g = geom;
			}

			Coordinate c = g.getEnvelopeInternal().centre();
			StringBuilder sb = new StringBuilder("https://maps.google.com/maps?ll=");
			sb.append(c.y);
			sb.append(",");
			sb.append(c.x);
			sb.append("&z=");
			sb.append(zoomlevel);
			sb.append("&t=");
			sb.append(layer);
			result = sb.toString();

			return result;

		}
		
		public static Geometry transform(Geometry geom, int srid) {
			if (geom == null) {
				return null;
			}

			Geometry result = GeometryTransform.transform(geom, srid);
			return result;
		}

		public static Geometry transform(Geometry geom, String to_proj) {
			if (geom == null) {
				return null;
			}

			Geometry result = GeometryTransform.transform(geom, to_proj);
			return result;
		}
		
	}
	
	
	
	
	
	//all functions that transform
	
	protected static Geom bind(Geometry geometry, int srid) {
		if (geometry == null) {
			return null;
		}
		if (srid == NO_SRID) {
			return new SimpleGeom(geometry);
		}
		return bind(geometry, SpatialReference.create(srid));
	}

	private static MapGeom bind(Geometry geometry, SpatialReference sr) {
		return new MapGeom(new MapGeometry(geometry, sr));
	}
	
	static class GeometryTransform {
		private static final CoordinateTransformFactory ctf = new CoordinateTransformFactory();
		private static final CRSFactory crsFactory = new CRSFactory();

		public GeometryTransform() {
		}

		public static Geometry transform(Geometry geom, int srid) {

			Geometry result = null;
			int geomsrid = geom.getSRID();

			CoordinateReferenceSystem fromcrs = crs(geomsrid);
			CoordinateReferenceSystem tocrs = crs(srid);
			CoordinateTransform trans = ctf.createTransform(fromcrs, tocrs);
			result = transformGeometry(trans, geom);

			return result;

		}

		public static Geometry transform(Geometry geom, String to_proj) {

			Geometry result = null;
			int geomsrid = geom.getSRID();

			CoordinateReferenceSystem fromcrs = crs(geomsrid);
			CoordinateReferenceSystem tocrs = crs(to_proj);
			CoordinateTransform trans = ctf.createTransform(fromcrs, tocrs);
			result = transformGeometry(trans, geom);

			return result;

		}

		protected static CoordinateReferenceSystem crs(int srid) {

			CoordinateReferenceSystem result = crsFactory.createFromName(String.format("epsg:%s", srid));

			return result;

		}

		protected static CoordinateReferenceSystem crs(String proj) {

			CoordinateReferenceSystem result = crsFactory.createFromParameters(null, proj);

			return result;

		}

		protected static Coordinate[] convert(ProjCoordinate[] projCoords) {
			Coordinate[] jtsCoords = new Coordinate[projCoords.length];
			for (int i = 0; i < projCoords.length; ++i) {
				jtsCoords[i] = new Coordinate(projCoords[i].x, projCoords[i].y);
			}
			return jtsCoords;
		}

		protected static ProjCoordinate[] convert(Coordinate[] jtsCoords) {
			ProjCoordinate[] projCoords = new ProjCoordinate[jtsCoords.length];
			for (int i = 0; i < jtsCoords.length; ++i) {
				projCoords[i] = new ProjCoordinate(jtsCoords[i].x, jtsCoords[i].y);
			}
			return projCoords;
		}

		protected static Coordinate[] transformCoordinates(CoordinateTransform ct, Coordinate[] in) {
			return convert(transformCoordinates(ct, convert(in)));
		}

		protected static ProjCoordinate[] transformCoordinates(CoordinateTransform ct, ProjCoordinate[] in) {
			ProjCoordinate[] out = new ProjCoordinate[in.length];
			for (int i = 0; i < in.length; ++i) {
				out[i] = ct.transform(in[i], new ProjCoordinate());
			}
			return out;
		}

		protected static Polygon transformPolygon(CoordinateTransform ct, Polygon polygon) {
			return polygon.getFactory().createPolygon(transformCoordinates(ct, polygon.getCoordinates()));
		}

		protected static Geometry transformPoint(CoordinateTransform ct, Point point) {
			return point.getFactory().createPoint(transformCoordinates(ct, point.getCoordinates())[0]);
		}

		protected static Geometry transformLinearRing(CoordinateTransform ct, LinearRing linearRing) {
			return linearRing.getFactory().createLinearRing(transformCoordinates(ct, linearRing.getCoordinates()));
		}

		protected static Geometry transformLineString(CoordinateTransform ct, LineString lineString) {
			return lineString.getFactory().createLineString(transformCoordinates(ct, lineString.getCoordinates()));
		}

		protected static Geometry transformMultiPolygon(CoordinateTransform ct, MultiPolygon multiPolygon) {
			Polygon[] polygon = new Polygon[multiPolygon.getNumGeometries()];
			for (int i = 0; i < polygon.length; ++i) {
				polygon[i] = multiPolygon.getFactory()
						.createPolygon(transformCoordinates(ct, multiPolygon.getGeometryN(i).getCoordinates()));
			}
			return multiPolygon.getFactory().createMultiPolygon(polygon);
		}

		protected static Geometry transformMultiPoint(CoordinateTransform ct, MultiPoint multiPoint) {
			return multiPoint.getFactory().createMultiPoint(transformCoordinates(ct, multiPoint.getCoordinates()));
		}

		protected static Geometry transformMultiLineString(CoordinateTransform ct, MultiLineString multiLineString) {
			LineString[] lineString = new LineString[multiLineString.getNumGeometries()];
			for (int i = 0; i < lineString.length; ++i) {
				lineString[i] = multiLineString.getFactory()
						.createLineString(transformCoordinates(ct, multiLineString.getGeometryN(i).getCoordinates()));
			}
			return multiLineString.getFactory().createMultiLineString(lineString);
		}

		protected static Geometry transformGeometryCollection(CoordinateTransform ct,
				GeometryCollection geometryCollection) {
			Geometry[] geometry = new Geometry[geometryCollection.getNumGeometries()];
			for (int i = 0; i < geometry.length; ++i) {
				geometry[i] = transformGeometry(ct, geometryCollection.getGeometryN(i));
			}
			return geometryCollection.getFactory().createGeometryCollection(geometry);
		}

		protected static Geometry transformGeometry(CoordinateTransform ct, Geometry geom) {
			if (geom instanceof Polygon) {
				return transformPolygon(ct, (Polygon) geom);
			} else if (geom instanceof Point) {
				return transformPoint(ct, (Point) geom);
			} else if (geom instanceof LinearRing) {
				return transformLinearRing(ct, (LinearRing) geom);
			} else if (geom instanceof LineString) {
				return transformLineString(ct, (LineString) geom);
			} else if (geom instanceof MultiPolygon) {
				return transformMultiPolygon(ct, (MultiPolygon) geom);
			} else if (geom instanceof MultiPoint) {
				return transformMultiPoint(ct, (MultiPoint) geom);
			} else if (geom instanceof MultiLineString) {
				return transformMultiLineString(ct, (MultiLineString) geom);
			} else if (geom instanceof GeometryCollection) {
				return transformGeometryCollection(ct, (GeometryCollection) geom);
			}
			return null;
		}
	}

	public static class SpatialReference {
		private static final CRSFactory crsFactory = new CRSFactory();
		public static CoordinateReferenceSystem crs;
		public static int wkid;

		public SpatialReference() {
		}

		public static SpatialReference create(int srid) {
			if (srid <= 0)
				throw new IllegalArgumentException("Invalid or unsupported wkid: " + wkid);

			SpatialReference sr = new SpatialReference();
			wkid = srid;
			crs = crsFactory.createFromName(String.format("epsg:%s", srid));

			return sr;
		}

		public int GetSRID() {
			return wkid;
		}
	}

	public enum GeometryType {
		// UNKNOWN(0),
		GEOMETRY(0), POINT(1), LINESTRING(2), POLYGON(3), MULTIPOINT(4), MULTILINESTRING(5), MULTIPOLYGON(
				6), GEOMETRYCOLLECTION(7);

		private static final Map<Integer, GeometryType> lookup = new HashMap<Integer, GeometryType>();

		static {
			for (GeometryType s : EnumSet.allOf(GeometryType.class))
				lookup.put(s.getCode(), s);
		}

		private int code;

		private GeometryType(int code) {
			this.code = code;
		}

		public int getCode() {
			return code;
		}

		public static GeometryType get(int code) {
			return lookup.get(code);
		}
	}

	// // GeometryType Enum - END
	public static interface Geom {
		Geometry g();

		SpatialReference sr();

		Geom wrap(Geometry g);
	}

	// interface Geom - END
	static class MapGeometry {
		Geometry m_geometry = null;
		SpatialReference sr = null;

		public MapGeometry(Geometry g, SpatialReference _sr) {
			m_geometry = g;

			sr = _sr;
			m_geometry.setSRID(getSRID());
		}

		public Geometry getGeometry() {
			return m_geometry;
		}

		public void setGeometry(Geometry geometry) {
			this.m_geometry = geometry;
		}

		public void setSpatialReference(SpatialReference _sr) {
			this.sr = _sr;
		}

		public SpatialReference getSpatialReference() {
			return sr;
		}

		public int getSRID() {
			return this.sr.GetSRID();
		}
	}

	static class MapGeom implements Geom {
		final MapGeometry mg;

		MapGeom(MapGeometry mg) {
			this.mg = Objects.requireNonNull(mg);
		}

		@Override
		public String toString() {
			return mg.toString();
		}

		public Geometry g() {
			return mg.getGeometry();
		}

		public Geom wrap(Geometry g) {
			return bind(g, this.mg.getSpatialReference());
		}

		public SpatialReference sr() {

			return mg.getSpatialReference();
		}

	}

	static class SimpleGeom implements Geom {
		final Geometry g;

		SimpleGeom(Geometry g) {
			this.g = Objects.requireNonNull(g);
		}

		@Override
		public String toString() {
			return g.toString();
		}

		public Geometry g() {
			return g;
		}

		public Geom wrap(Geometry g) {
			return new SimpleGeom(g);
		}

		public SpatialReference sr() {
			// TODO Auto-generated method stub
			return SPATIAL_REFERENCE;
		}

	}
	

}
