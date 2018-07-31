                  MOUNT WILSON SUNSPOT DATA NOTES

Sunspot umbral position and area information were digitized from
the Mount Wilson daily white-light solar images some years ago
(Howard, Gilman, and Gilman, 1984). These photographic images
exist in a series that extends from 1917 through the present
time. The digitized data extend from 1917 through 1985. Details
about the observations, the measurement procedure, and the
analysis techniques used earlier may be obtained from the earlier
reference (Howard, Gilman, and Gilman, 1984).

I. RAW DATA FILES

These files contain the raw data measurements made in the early
1980's. The digital measurements were made manually with a
digitizing pad and cursor. The image was aligned so that the axis
of rotation was parallel to the Y axis of the measuring pad. Each
plate was aligned for the position measurements at approximately
the same location on the measuring pad. For each day's data
there are eight equally-spaced limb measurements in X and Y. For
each sunspot there are two position measurements (in X and Y).
From these two measurements it is possible to determine the
sunspot umbral area (Howard, Gilman, and Gilman, 1984). On the
measuring pad, the arbitrary length units are lower in the
(heliocentric) south and lower in the (heliocentric) west. All
units are positive.

Various other information is included in the records for each
day. The details for these are listed below.

The user should be aware that these data do not contain any
corrections for atmospheric refraction (Howard, Gilman, and
Gilman, 1984) or any other effects, such as the optical
aberrations corrected more recently (Howard, Gupta, and
Sivaraman, 1999). These are pure raw data files intended for
those users who wish to make their own corrections.

The records are formatted as follows:

    format 8I2,I1,I3
 I2     Day of the month
 I2     Month
 I2     Year (2-digit, 17 through 85)
 I2     Hour of observation (This is Pacific Standard Time, which
         is 8 hours west of Greenwich)
 I2     Minute of observation
 I2     Estimate of quality of the image (0 to 10; high numbers
         represent better quality)
 I2     Estimate of average photographic image density on the plate
 I2     Estimate of photographic density off the image
 I1     Obscuration (0 = none; 1 = image partly obscured)
 I3     NSPOT (Number of sunspots for this day's observation)

    format 2F4.1,F5.1,3F7.3/8F7.3/8F7.3
 F4.1   Estimate of seeing from the image (Low numbers are poor
         seeing)
 F4.1   Temperature in the telescope observing room (unheated)
         (In Celsius before 1921; in Fahrenheit from 1921 on)
 F5.1   Telescope focus setting (arbitrary)
 F7.1   X center position from first rough circle solution (no
         corrections applied to positions)
 F7.1   Y center position from first rough circle solution (no
         corrections applied to positions)
 F7.2   Radius determination from first rough circle solution (no
         corrections applied to positions)
 8F7.3  X limb positions (no corrections applied to positions)
 8f7.3  Y limb positions (corresponding in location in the array
         to the X limb positions above) (no corrections applied
         to positions)

    format 4A1
 4A1    Initials of observer

    format 10F7.3 (This record is present only if NSPOT is
            greater than 0)
 10F7.3    Array of X points of sunspot umbral measurements
            (there are NSPOT of these)
 10F7.3    Array of Y points of sunspot umbral measurements
            (there are NSPOT of these, corresponding in location
            in the array to the X sunspot positions above)

II. SUNSPOT GROUP DATA

These files contain reduced data for sunspot groups as defined in
the earlier publications (Howard, Gilman, and Gilman 1984;
Howard, Gupta, and Sivaraman 1999). In general, as is pointed out
in these publications, we were very conservative in our selection
of groups (but especially for sunspot returns), preferring to
miss some real next-day returns rather than take a chance on
including spurious data.

The data for each group comes from observations from two
consecutive days. Thus we can include rotation and meridional
drift results and other quantities that may vary with time. There
is no lifetime information in this data set because we did not
keep track of any group for more than two consecutive days.

The reader should refer to the papers cited below for details
concerning the various quantities listed here. Latitudes and
longitudes are heliocentric. Negative latitudes are in the south;
negative longitudes are in the east. Areas are in micro-
hemispheres. For definitions of 'leading' and 'following'
sunspots, see Howard, Gilman, and Gilman, 1984.

The records are formatted as follows:

    format 9I2,2F10.4
 I2     Month
 I2     Day of the month
 I2     Year (2-digit, 17 through 85)
 I2     Hour of observation (This is Pacific Standard Time, which
         is 8 hours west of Greenwich)
 I2     Minute of observation
 I2     Number of sunspots in the leading portion of the group
 I2     Number of sunspots in the following portion of the group
 I2     Total number of sunspots in the group on day 1
 I2     Total number of sunspots in the group on day 2
 F10.4  Day number of first observation
 F10.4  Day number of second observation
        (These day numbers are in a system that starts Jan 1, 1915)

    format 9F7.3/8F7.3,F7.2
 F7.3   Rotation rate of the group in deg/day sidereal
 F7.3   Latitude drift of the group in deg/day (+ to north)
 F7.3   Area-weighted latitude of the group on day 1
 F7.3   Area-weighted longitude of the group on day 1
 F7.3   Area-weighted longitude of the group on day 2
 F7.3   Area of the group on day 1
 F7.3   Area of the group on day 2
 F7.3   Tilt angle of the group on day 1 in degrees (1)
 F7.3   Change (day2 - day1) in the polarity separation (1)
 F7.3   Area-weighted longitude of the sunspots in the leading portion (1)
 F7.3   Area-weighted longitude of the sunspots in the following portion (1)
 F7.3   Area-weighted latitude of the sunspots in the leading portion (1)
 F7.3   Area-weighted latitude of the sunspots in the following portion (1)
 F7.3   Total area of the leading sunspots (1)
 F7.3   Total area of the following sunspots (1)
 F7.3   Area-weighted longitude of leading sunspots on day 2 (2)
 F7.3   Area-weighted longitude of following sunspots on day 2 (2)
 F7.2   Difference in the tilt angles (day2 - day1) in deg/day (1) (2)

(1) This quantity is zero when there is only one spot in the
group on day 1.
(2) This quantity is zero when there is only one spot in the
group on day 2.

III. INDIVIDUAL SUNSPOT DATA

These data refer to individual sunspots in the data set. Again
there is no lifetime information. Also, as with the groups, the
selection criteria for identifying a sunspot as a 'return' on the
next day have been made very conservative (see the papers
referenced below), so that it is likely that some sunspots were
missed in this selection, but this reduces the chances of
mis-identifications.

The records are formatted as follows:

    format 9I2,2F10.4
 I2     Month
 I2     Day of the month
 I2     Year (2-digit, 17 through 85)
 I2     Hour of observation (This is Pacific Standard Time, which
         is 8 hours west of Greenwich)
 I2     Minute of observation
 I2     Number of leading sunspots in the group on day 1
 I2     Number of following sunspots in the group on day 1
 I2     Total number of sunspots in the group on day 1
 I2     Total number of sunspots in the group on day 2
 F10.4  Day number of first observation
 F10.4  Day number of second observation
        (These day numbers are in a system that starts Jan 1, 1915)

    format 9F7.3/8F7.3,F7.2
 F7.3   Rotation rate of the sunspot in deg/day sidereal
 F7.3   Latitude drift of the sunspot in deg/day (+ to north)
 F7.3   Latitude of the sunspot on day 1
 F7.3   Latitude of the sunspot on day 2
 F7.3   Longitude of the sunspot on day 1
 F7.3   longitude of the sunspot on day 2
 F7.3   Area of the sunspot on day 1
 F7.3   Area of the sunspot on day 2
 F7.3   Latitude of the sunspot group on day 1
 F7.3   Longitude of the sunspot group on day 1
 F7.3   Area of the sunspot group on day 1
 F7.3   Area of the sunspot group on day 2
 F7.3   Polarity separation of the group on day 1
 F7.3   Polarity separation change of the group (day 2 - day 1)
 F7.3   Tilt angle of the group on day 1
 F7.3   Tilt angle of the group on day 2
 F7.3   Rotation rate of the group in deg/day sidereal
 F7.3   Latitude drift of the group in deg/day

References

Howard, Gilman, and Gilman, 1984, ApJ, 283, 373-384.

Howard, Gupta, and Sivaraman, 1999, Solar Physics, 186, 25-41.
