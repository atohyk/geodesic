#include <stdio.h>
#include <math.h>
#include <time.h>

const double earth_radius = 6378137.0;//in metres

double rad2deg(double rad){
	return (rad * 180.0 / M_PI);
}

double deg2rad(double deg){
	return (deg * M_PI / 180.0);
}

//Uses simple pythagoras with small angle approximation.
//Should only be used closed to the equator for small distances.
double geodesic_pythagoras(double lat1, double long1, double lat2, double long2){
	return sqrt(pow((long2 - long1) * earth_radius, 2) + pow((lat2 - lat1) * earth_radius, 2));
}

//most accurate, should be used for long distances
//takes 500us on a Cortex M4
double geodesic_vincenty(double lat1, double long1, double lat2, double long2){
	//init
	const double a = earth_radius;
	const double f = 1/298.257223563;
	const double b = 6356752.314245;
	const double threshold = 10E-12;
	double tanu1 = (1-f)*tan(lat1);
	double cosu1 = 1/sqrt(1+pow(tanu1, 2));
	double sinu1 = tanu1*cosu1;
	double tanu2 = (1-f)*tan(lat2);
	double cosu2 = 1/sqrt(1+pow(tanu2, 2));
	double sinu2 = tanu2*cosu2;
	double L = long2-long1;
	double lambda = L;

	//loop vars
	double sigma = 0, C, sinSigma = 0, cosSigma = 0, sinAlpha, sinLambda;
	double cosLambda, cos2Alpha = 0, cos2SigM = 0, oldLambda = lambda+1;

	while(fabs(oldLambda - lambda) > threshold){
		oldLambda = lambda;
		sinLambda = sin(lambda);
		cosLambda = cos(lambda);
		sinSigma = sqrt(pow(cosu2*sinLambda, 2) + pow(cosu1*sinu2-sinu1*cosu2*cosLambda, 2));
		cosSigma = sinu1*sinu2+cosu1*cosu2*cosLambda;
		sigma = atan2(sinSigma, cosSigma);
		sinAlpha = cosu1*cosu2*sinLambda/sinSigma;
		cos2Alpha = 1-pow(sinAlpha, 2);
		cos2SigM = cosSigma-2*sinu1*sinu2/cos2Alpha;
		C = f/16*cos2Alpha*(4+f*(4-3*cos2Alpha));
		lambda = L+(1-C)*f*sinAlpha*(sigma+C*sinAlpha*(cos2SigM+C*cosSigma*(-1+2*cos2SigM)));
	}

	double usq = cos2Alpha * (pow(a,2) - pow(b,2))/pow(b,2);
	double A = 1 + usq/16384 * (4096 + usq * (-768 + usq*(320-175*usq)));
	double B = usq/1024*(256+usq*(-128+usq*(74-47*usq)));
	double dSigma = B*sinSigma*(cos2SigM+B/4*(cosSigma*(-1+2*pow(cos2SigM,2))-
			B/6*cos2SigM*(-3+4*pow(sinSigma, 2))*(-3+4*pow(cos2SigM, 2))));
	double s = b*A*(sigma-dSigma);
	return s;
}

//assumes great circle 
double geodesic_haversine(double lat1, double long1, double lat2, double long2){
	double deltaLong = long2 - long1;
	double deltaSigma = atan2(sqrt(pow(cos(lat2)*sin(deltaLong), 2)+
		pow(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(deltaLong), 2)),
		sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(deltaLong));
	return earth_radius*deltaSigma;
}

double geodesic_cos_long(double lat1, double long1, double lat2, double long2){
	return sqrt(pow((cos((lat1+lat2)/2)*(long2-long1)) * earth_radius, 2) + pow((lat2 - lat1) * earth_radius, 2));
}

int main(){
	//in degrees
	double deglat1 = 1.293642, deglong1 = 103.857045; //suntec city
	double deglat2 = 1.304369, deglong2 = 103.831974; //ion orchard
	//should be around 3.03 km

	double lat1 = deg2rad(deglat1);
	double lat2 = deg2rad(deglat2);
	double long1 = deg2rad(deglong1);
	double long2 = deg2rad(deglong2);

	printf("=====TEST CASE 1=====\n");
	printf("Lat 1: %f  Long 1: %f\nLat 2: %f  Long2: %f\n", lat1, long1, lat2, long2);
	printf("Naive: %f\n", geodesic_pythagoras(lat1, long1, lat2, long2));
	printf("CosLong: %f\n", geodesic_cos_long(lat1, long1, lat2, long2));
	printf("Haversine: %f\n", geodesic_haversine(lat1, long1, lat2, long2));
	printf("Vincenty: %f\n", geodesic_vincenty(lat1, long1, lat2, long2));
	printf("=====END TEST CASE=====\n\n");

	//Changi Airport
	deglat1 = 1.360307;
	deglong1 = 103.989809;
	//Taiwan Taoyuan International Airport
	deglat2 = 25.076984;
	deglong2 = 121.230774;

	lat1 = deg2rad(deglat1);
	lat2 = deg2rad(deglat2);
	long1 = deg2rad(deglong1);
	long2 = deg2rad(deglong2);

	printf("=====TEST CASE 2=====\n");
	printf("Lat 1: %f  Long 1: %f\nLat 2: %f  Long2: %f\n", lat1, long1, lat2, long2);
	printf("Naive: %f\n", geodesic_pythagoras(lat1, long1, lat2, long2));
	printf("CosLong: %f\n", geodesic_cos_long(lat1, long1, lat2, long2));
	printf("Haversine: %f\n", geodesic_haversine(lat1, long1, lat2, long2));
	printf("Vincenty: %f\n", geodesic_vincenty(lat1, long1, lat2, long2));
	printf("=====END TEST CASE=====\n\n");
}