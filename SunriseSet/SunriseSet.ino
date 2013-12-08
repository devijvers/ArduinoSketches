
#include <stdio.h>
#include <math.h>
#include <Time.h>
#include <Wire.h>  
#include <DS1307RTC.h>  // a basic DS1307 library that returns time as a time_t


#define RELAY_PIN 12
// see http://itouchmap.com/latlong.html to find lat/long for your position
#define HOME_LONGITUDE 4.787389
#define HOME_LATITUDE  50.850214

// This algorithm only support DST changes on a Sunday
// Sorry Israel, Iran, Palestina and Syria.
#define DST_ENABLED true
#define DST_




/* This macro computes times for sunrise/sunset.                      */
/* Sunrise/set is considered to occur when the Sun's upper limb is    */
/* 35 arc minutes below the horizon (this accounts for the refraction */
/* of the Earth's atmosphere).                                        */
#define sun_rise_set(year,month,day,lon,lat,rise,set)  \
        __sunriset__( year, month, day, lon, lat, -35.0/60.0, 1, rise, set )



int weekdayOffsets[] = {-1, 0, 6, 5, 4, 3, 2, 1};

int __getDSTDay__(int _month, int _year) {
  tmElements_t tm;
  tm.Second = 0;
  tm.Minute = 0;
  tm.Hour = 0;
  tm.Day = 25;
  tm.Month = _month;
  tm.Year = _year - 1970;

  time_t d = makeTime(tm);
  int _weekday = weekday(d);
  if (_weekday == 1) return 25;
  int offset = weekdayOffsets[_weekday];
  return offset + 25;
}

#define getMarchDSTStartDay(_year) \
    __getDSTDay__(3, _year)

#define getOctoberDSTEndDay(_year) \
    __getDSTDay__(10, _year)

boolean isDaylightSavingTime(time_t t) {
  if (DST_ENABLED) {
    int _dayOfMonth = day(t);
    int _month = month(t);
    int _year = year(t);
    int dstDay;
    
    switch (_month) {
       case 1: 
       case 2: 
         return false;
       case 3:
         dstDay = getMarchDSTStartDay(_year);
         return (_dayOfMonth >= dstDay);
       case 4:
       case 5:
       case 6:
       case 7:
       case 8:
       case 9: 
         return true;
       case 10:
         dstDay = getOctoberDSTEndDay(_year);
         return (_dayOfMonth < dstDay);
    }
    
    return false;
  } else {
    return false;
  }
}

void setup() {
  Serial.begin(9600);
  pinMode(RELAY_PIN, OUTPUT);
  pinMode(13, OUTPUT);
  
  digitalWrite(RELAY_PIN, HIGH);
  setSyncProvider(RTC.get);   // the function to get the time from the RTC
  if(timeStatus()!= timeSet) {
      digitalWrite(13, LOW);
  } else {
      digitalWrite(13, HIGH);
  }
}

void loop() {

  tmElements_t tm;
  tm.Second = 0;
  tm.Minute = 0;
  tm.Hour = 0;
  tm.Day = 6;
  tm.Month = 9;
  tm.Year = 2013 - 1970;

  
  time_t t = RTC.get();
  boolean dst = isDaylightSavingTime(t);
  int _minute = minute(t);
  int _hour = hour(t);
  int _year = year(t);
  int _month = month(t);
  int _day = day(t);
  double _rise, _set;
  
  sun_rise_set(_year,_month,_day,HOME_LONGITUDE,HOME_LATITUDE,&_rise,&_set);
  _rise += dst ? 2.0 : 1.0;
  _set += dst ? 2.0 : 1.0;
  
  double n = double(_hour) + double(_minute) / 60;
  Serial.print("Now: ");
  Serial.print(n);
  Serial.print(", Sunrise: ");
  Serial.print(_rise);
  Serial.print(", sunset: ");
  Serial.println(_set);
  
  if (n >= _set || (n >= 0.0 && n <= 1.0) || (n >= 5.0 && n <= _rise)) {
    digitalWrite(RELAY_PIN, LOW);
  } else {
    digitalWrite(RELAY_PIN, HIGH);
  }

  delay(2000);   
}



/* +++Date last modified: 05-Jul-1997 */
/* Updated comments, 05-Aug-2013 */

/*

SUNRISET.C - computes Sun rise/set times, start/end of twilight, and
             the length of the day at any date and latitude

Written as DAYLEN.C, 1989-08-16

Modified to SUNRISET.C, 1992-12-01

(c) Paul Schlyter, 1989, 1992

Released to the public domain by Paul Schlyter, December 1992

*/


/* A macro to compute the number of days elapsed since 2000 Jan 0.0 */
/* (which is equal to 1999 Dec 31, 0h UT)                           */

#define days_since_2000_Jan_0(y,m,d) \
    (367L*(y)-((7*((y)+(((m)+9)/12)))/4)+((275*(m))/9)+(d)-730530L)

/* Some conversion factors between radians and degrees */

#ifndef PI
 #define PI        3.1415926535897932384
#endif

#define RADEG     ( 180.0 / PI )
#define DEGRAD    ( PI / 180.0 )

/* The trigonometric functions in degrees */

#define sind(x)  sin((x)*DEGRAD)
#define cosd(x)  cos((x)*DEGRAD)
#define tand(x)  tan((x)*DEGRAD)

#define atand(x)    (RADEG*atan(x))
#define asind(x)    (RADEG*asin(x))
#define acosd(x)    (RADEG*acos(x))
#define atan2d(y,x) (RADEG*atan2(y,x))

/* Function prototypes */

double __daylen__( int _year, int _month, int _day, double lon, double lat,
                   double altit, int upper_limb );

int __sunriset__( int _year, int _month, int _day, double lon, double lat,
                  double altit, int upper_limb, double *rise, double *set );



void sunpos( double d, double *lon, double *r );

void sun_RA_dec( double d, double *RA, double *dec, double *r );

double revolution( double x );

double rev180( double x );

double GMST0( double d );

/* The "workhorse" function for sun rise/set times */

int __sunriset__( int _year, int _month, int _day, double lon, double lat,
                  double altit, int upper_limb, double *trise, double *tset )
/***************************************************************************/
/* Note: year,month,date = calendar date, 1801-2099 only.             */
/*       Eastern longitude positive, Western longitude negative       */
/*       Northern latitude positive, Southern latitude negative       */
/*       The longitude value IS critical in this function!            */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight.                   */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing rise/set     */
/*               times, and to zero when computing start/end of       */
/*               twilight.                                            */
/*        *rise = where to store the rise time                        */
/*        *set  = where to store the set  time                        */
/*                Both times are relative to the specified altitude,  */
/*                and thus this function can be used to compute       */
/*                various twilight times, as well as rise/set times   */
/* Return value:  0 = sun rises/sets this day, times stored at        */
/*                    *trise and *tset.                               */
/*               +1 = sun above the specified "horizon" 24 hours.     */
/*                    *trise set to time when the sun is at south,    */
/*                    minus 12 hours while *tset is set to the south  */
/*                    time plus 12 hours. "Day" length = 24 hours     */
/*               -1 = sun is below the specified "horizon" 24 hours   */
/*                    "Day" length = 0 hours, *trise and *tset are    */
/*                    both set to the time when the sun is at south.  */
/*                                                                    */
/**********************************************************************/
{
      double  d,  /* Days since 2000 Jan 0.0 (negative before) */
      sr,         /* Solar distance, astronomical units */
      sRA,        /* Sun's Right Ascension */
      sdec,       /* Sun's declination */
      sradius,    /* Sun's apparent radius */
      t,          /* Diurnal arc */
      tsouth,     /* Time when Sun is at south */
      sidtime;    /* Local sidereal time */

      int rc = 0; /* Return cde from function - usually 0 */

      /* Compute d of 12h local mean solar time */
      d = days_since_2000_Jan_0(_year,_month,_day) + 0.5 - lon/360.0;

      /* Compute the local sidereal time of this moment */
      sidtime = revolution( GMST0(d) + 180.0 + lon );

      /* Compute Sun's RA, Decl and distance at this moment */
      sun_RA_dec( d, &sRA, &sdec, &sr );

      /* Compute time when Sun is at south - in hours UT */
      tsouth = 12.0 - rev180(sidtime - sRA)/15.0;

      /* Compute the Sun's apparent radius in degrees */
      sradius = 0.2666 / sr;

      /* Do correction to upper limb, if necessary */
      if ( upper_limb )
            altit -= sradius;

      /* Compute the diurnal arc that the Sun traverses to reach */
      /* the specified altitude altit: */
      {
            double cost;
            cost = ( sind(altit) - sind(lat) * sind(sdec) ) /
                  ( cosd(lat) * cosd(sdec) );
            if ( cost >= 1.0 )
                  rc = -1, t = 0.0;       /* Sun always below altit */
            else if ( cost <= -1.0 )
                  rc = +1, t = 12.0;      /* Sun always above altit */
            else
                  t = acosd(cost)/15.0;   /* The diurnal arc, hours */
      }

      /* Store rise and set times - in hours UT */
      *trise = tsouth - t;
      *tset  = tsouth + t;

      return rc;
}  /* __sunriset__ */



/* The "workhorse" function */


double __daylen__( int _year, int _month, int _day, double lon, double lat,
                   double altit, int upper_limb )
/**********************************************************************/
/* Note: year,month,date = calendar date, 1801-2099 only.             */
/*       Eastern longitude positive, Western longitude negative       */
/*       Northern latitude positive, Southern latitude negative       */
/*       The longitude value is not critical. Set it to the correct   */
/*       longitude if you're picky, otherwise set to to, say, 0.0     */
/*       The latitude however IS critical - be sure to get it correct */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight.                   */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing day length   */
/*               and to zero when computing day+twilight length.      */
/**********************************************************************/
{
      double  d,  /* Days since 2000 Jan 0.0 (negative before) */
      obl_ecl,    /* Obliquity (inclination) of Earth's axis */
      sr,         /* Solar distance, astronomical units */
      slon,       /* True solar longitude */
      sin_sdecl,  /* Sine of Sun's declination */
      cos_sdecl,  /* Cosine of Sun's declination */
      sradius,    /* Sun's apparent radius */
      t;          /* Diurnal arc */

      /* Compute d of 12h local mean solar time */
      d = days_since_2000_Jan_0(_year,_month,_day) + 0.5 - lon/360.0;

      /* Compute obliquity of ecliptic (inclination of Earth's axis) */
      obl_ecl = 23.4393 - 3.563E-7 * d;

      /* Compute Sun's ecliptic longitude and distance */
      sunpos( d, &slon, &sr );

      /* Compute sine and cosine of Sun's declination */
      sin_sdecl = sind(obl_ecl) * sind(slon);
      cos_sdecl = sqrt( 1.0 - sin_sdecl * sin_sdecl );

      /* Compute the Sun's apparent radius, degrees */
      sradius = 0.2666 / sr;

      /* Do correction to upper limb, if necessary */
      if ( upper_limb )
            altit -= sradius;

      /* Compute the diurnal arc that the Sun traverses to reach */
      /* the specified altitude altit: */
      {
            double cost;
            cost = ( sind(altit) - sind(lat) * sin_sdecl ) /
                  ( cosd(lat) * cos_sdecl );
            if ( cost >= 1.0 )
                  t = 0.0;                      /* Sun always below altit */
            else if ( cost <= -1.0 )
                  t = 24.0;                     /* Sun always above altit */
            else  t = (2.0/15.0) * acosd(cost); /* The diurnal arc, hours */
      }
      return t;
}  /* __daylen__ */


/* This function computes the Sun's position at any instant */

void sunpos( double d, double *lon, double *r )
/******************************************************/
/* Computes the Sun's ecliptic longitude and distance */
/* at an instant given in d, number of days since     */
/* 2000 Jan 0.0.  The Sun's ecliptic latitude is not  */
/* computed, since it's always very near 0.           */
/******************************************************/
{
      double M,         /* Mean anomaly of the Sun */
             w,         /* Mean longitude of perihelion */
                        /* Note: Sun's mean longitude = M + w */
             e,         /* Eccentricity of Earth's orbit */
             E,         /* Eccentric anomaly */
             x, y,      /* x, y coordinates in orbit */
             v;         /* True anomaly */

      /* Compute mean elements */
      M = revolution( 356.0470 + 0.9856002585 * d );
      w = 282.9404 + 4.70935E-5 * d;
      e = 0.016709 - 1.151E-9 * d;

      /* Compute true longitude and radius vector */
      E = M + e * RADEG * sind(M) * ( 1.0 + e * cosd(M) );
            x = cosd(E) - e;
      y = sqrt( 1.0 - e*e ) * sind(E);
      *r = sqrt( x*x + y*y );              /* Solar distance */
      v = atan2d( y, x );                  /* True anomaly */
      *lon = v + w;                        /* True solar longitude */
      if ( *lon >= 360.0 )
            *lon -= 360.0;                   /* Make it 0..360 degrees */
}

void sun_RA_dec( double d, double *RA, double *dec, double *r )
/******************************************************/
/* Computes the Sun's equatorial coordinates RA, Decl */
/* and also its distance, at an instant given in d,   */
/* the number of days since 2000 Jan 0.0.             */
/******************************************************/
{
      double lon, obl_ecl, x, y, z;

      /* Compute Sun's ecliptical coordinates */
      sunpos( d, &lon, r );

      /* Compute ecliptic rectangular coordinates (z=0) */
      x = *r * cosd(lon);
      y = *r * sind(lon);

      /* Compute obliquity of ecliptic (inclination of Earth's axis) */
      obl_ecl = 23.4393 - 3.563E-7 * d;

      /* Convert to equatorial rectangular coordinates - x is unchanged */
      z = y * sind(obl_ecl);
      y = y * cosd(obl_ecl);

      /* Convert to spherical coordinates */
      *RA = atan2d( y, x );
      *dec = atan2d( z, sqrt(x*x + y*y) );

}  /* sun_RA_dec */


/******************************************************************/
/* This function reduces any angle to within the first revolution */
/* by subtracting or adding even multiples of 360.0 until the     */
/* result is >= 0.0 and < 360.0                                   */
/******************************************************************/

#define INV360    ( 1.0 / 360.0 )

double revolution( double x )
/*****************************************/
/* Reduce angle to within 0..360 degrees */
/*****************************************/
{
      return( x - 360.0 * floor( x * INV360 ) );
}  /* revolution */

double rev180( double x )
/*********************************************/
/* Reduce angle to within +180..+180 degrees */
/*********************************************/
{
      return( x - 360.0 * floor( x * INV360 + 0.5 ) );
}  /* revolution */


/*******************************************************************/
/* This function computes GMST0, the Greenwich Mean Sidereal Time  */
/* at 0h UT (i.e. the sidereal time at the Greenwhich meridian at  */
/* 0h UT).  GMST is then the sidereal time at Greenwich at any     */
/* time of the day.  I've generalized GMST0 as well, and define it */
/* as:  GMST0 = GMST - UT  --  this allows GMST0 to be computed at */
/* other times than 0h UT as well.  While this sounds somewhat     */
/* contradictory, it is very practical:  instead of computing      */
/* GMST like:                                                      */
/*                                                                 */
/*  GMST = (GMST0) + UT * (366.2422/365.2422)                      */
/*                                                                 */
/* where (GMST0) is the GMST last time UT was 0 hours, one simply  */
/* computes:                                                       */
/*                                                                 */
/*  GMST = GMST0 + UT                                              */
/*                                                                 */
/* where GMST0 is the GMST "at 0h UT" but at the current moment!   */
/* Defined in this way, GMST0 will increase with about 4 min a     */
/* day.  It also happens that GMST0 (in degrees, 1 hr = 15 degr)   */
/* is equal to the Sun's mean longitude plus/minus 180 degrees!    */
/* (if we neglect aberration, which amounts to 20 seconds of arc   */
/* or 1.33 seconds of time)                                        */
/*                                                                 */
/*******************************************************************/

double GMST0( double d )
{
      double sidtim0;
      /* Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  */
      /* L = M + w, as defined in sunpos().  Since I'm too lazy to */
      /* add these numbers, I'll let the C compiler do it for me.  */
      /* Any decent C compiler will add the constants at compile   */
      /* time, imposing no runtime or code overhead.               */
      sidtim0 = revolution( ( 180.0 + 356.0470 + 282.9404 ) +
                          ( 0.9856002585 + 4.70935E-5 ) * d );
      return sidtim0;
}  /* GMST0 */

