function [YY MO DD HH MI] = getSeasonStartDate( SeasonStr, YearIn )
%GETSEASONSTARTDATE load the year, month, day, hour, and minute of the
% start of that season in UT during 1900-2100
%
%
%   SYNTAX:
%                        getSeasonStartDate( SeasonStr, YearIn );
%     [YY MO DD HH MI] = getSeasonStartDate(...);
%
%   INPUT:
%     SeasonStr -   A character array specifying the season, either
%                   'spring'|'summer'|'fall'|'winter'.
%     YearIn    -   The year of the season change.
%
%   OUTPUT (optional):
%     YY        -   The year of the season change.
%     MO        -   The month of the season change.
%     DD        -   The day of the season change.
%     HH        -   The hour of the season change.
%     MI        -   The minute of the season change.
%
%   DESCRIPTION:
%     To find the the start of the requested season, which changes slightly
%     year to year.  The information for years 1900-2100 is taken from the
%     NASA Goddard Institute for Space Studies' Climate Model Simulations
%     <a href="http://aom.giss.nasa.gov/srvernal.html">site</a>
%     on 2011/08/31.
%
%     This function was first written to sort measurements by date to allow
%     me to accurately search for seasonal patterns.
%
%   SEE ALSO:
%     DATE, DATESTR, DATENUM, NOW
%     
%   ---
%   MFILE:   getSeasonStartDate.m
%   VERSION: 1.0 (2011/08/31) 
%   MATLAB:  7.8.0.347 (R2009a)
%   AUTHOR:  Anthony Bathgate
%   CONTACT: tony.bathgate@usask.ca
%
%   ADDITIONAL NOTES:
%     - In the future I might be persuaded to re-write this as a function
%     that calculates the dates on the fly, rather than using a lookup
%     table (because it would work for more years).
%
%   REVISIONS:
%   1.0      Released. (2011/8/31)
%   
%   DISCLAIMER:
%   getSeasonStartDate.m is provided "as is" without warranty of any kind, 
%   under the revised BSD license.
%
%   Copyright (C) 2011 Tony Bathgate, University of Saskatchewan, ISAS.



    % Verify the type
    if( ~ischar(SeasonStr) )
        error( ['first argument must be string ' ...
            '''spring''|''summer''|''fall''|''winter'''])
    end

    [dates] = LoadEquinoxSolsticeDates();

    % Ensure the year is in the table
    if( YearIn<dates(1,1) || YearIn>dates(end,1) )
        error( [ 'second argument must be a year in the set [' ...
            num2str( dates(1,1) ) ' ' num2str( dates(end,1) ) ']' ] );
    end

    % Find the row corresponding to the year
    yearIndex = YearIn-dates(1,1)+1;

    % Find the column corresponding to the season
    SeasonStr  = lower( SeasonStr );
    switch( SeasonStr )
        case 'spring'
            seasonIndex = 2;
        case 'summer'
            seasonIndex = 6;
        case 'fall'
            seasonIndex = 10;
        case 'winter'
            seasonIndex = 14;
        otherwise
            error(['first argument must be string ' ...
                '''spring''|''summer''|''fall''|''winter''']);
    end

    % Get the actual values
    YY  = dates(yearIndex,1);
    MO  = dates(yearIndex,seasonIndex+0);
    DD  = dates(yearIndex,seasonIndex+1);
    HH  = dates(yearIndex,seasonIndex+2);
    MI  = dates(yearIndex,seasonIndex+3);

end


function [dates] = LoadEquinoxSolsticeDates()
    % The Equinox and Solstice data is all taken from 
    % NASA Goddard Institute for Space Studies' Climate Model Simulations page
    % a href="http://aom.giss.nasa.gov/srvernal.html
    dates = ...
    [1900     3 21  1 30     6 21 21 30     9 23 12 04    12 22  6 32; ...
     1901     3 21  7 19     6 22  3 18     9 23 17 53    12 22 12 22; ...
     1902     3 21 13 08     6 22  9 06     9 23 23 42    12 22 18 12; ...
     1903     3 21 18 58     6 22 14 54     9 24  5 31    12 23  0 01; ...
     1904     3 21  0 47     6 21 20 43     9 23 11 19    12 22  5 51; ...
     1905     3 21  6 36     6 22  2 31     9 23 17 08    12 22 11 41; ...
     1906     3 21 12 25     6 22  8 19     9 23 22 57    12 22 17 31; ...
     1907     3 21 18 14     6 22 14 07     9 24  4 46    12 22 23 20; ...
     1908     3 21  0 04     6 21 19 55     9 23 10 34    12 22  5 10; ...
     1909     3 21  5 53     6 22  1 43     9 23 16 23    12 22 11 00; ...
     1910     3 21 11 42     6 22  7 31     9 23 22 12    12 22 16 50; ...
     1911     3 21 17 31     6 22 13 20     9 24  4 01    12 22 22 39; ...
     1912     3 20 23 20     6 21 19 08     9 23  9 49    12 22  4 29; ...
     1913     3 21  5 10     6 22  0 56     9 23 15 38    12 22 10 19; ...
     1914     3 21 10 59     6 22  6 44     9 23 21 27    12 22 16 09; ...
     1915     3 21 16 48     6 22 12 32     9 24  3 15    12 22 21 59; ...
     1916     3 20 22 37     6 21 18 20     9 23  9 04    12 22  3 48; ...
     1917     3 21  4 26     6 22  0 08     9 23 14 53    12 22  9 38; ...
     1918     3 21 10 16     6 22  5 57     9 23 20 42    12 22 15 28; ...
     1919     3 21 16 05     6 22 11 45     9 24  2 30    12 22 21 18; ...
     1920     3 20 21 54     6 21 17 33     9 23  8 19    12 22  3 07; ...
     1921     3 21  3 43     6 21 23 21     9 23 14 08    12 22  8 57; ...
     1922     3 21  9 32     6 22  5 09     9 23 19 57    12 22 14 47; ...
     1923     3 21 15 22     6 22 10 57     9 24  1 45    12 22 20 37; ...
     1924     3 20 21 11     6 21 16 45     9 23  7 34    12 22  2 26; ...
     1925     3 21  3 00     6 21 22 34     9 23 13 23    12 22  8 16; ...
     1926     3 21  8 49     6 22  4 22     9 23 19 12    12 22 14 06; ...
     1927     3 21 14 38     6 22 10 10     9 24  1 00    12 22 19 56; ...
     1928     3 20 20 28     6 21 15 58     9 23  6 49    12 22  1 45; ...
     1929     3 21  2 17     6 21 21 46     9 23 12 38    12 22  7 35; ...
     1930     3 21  8 06     6 22  3 34     9 23 18 26    12 22 13 25; ...
     1931     3 21 13 55     6 22  9 22     9 24  0 15    12 22 19 15; ...
     1932     3 20 19 44     6 21 15 11     9 23  6 04    12 22  1 04; ...
     1933     3 21  1 34     6 21 20 59     9 23 11 53    12 22  6 54; ...
     1934     3 21  7 23     6 22  2 47     9 23 17 41    12 22 12 44; ...
     1935     3 21 13 12     6 22  8 35     9 23 23 30    12 22 18 34; ...
     1936     3 20 19 01     6 21 14 23     9 23  5 19    12 22  0 23; ...
     1937     3 21  0 50     6 21 20 11     9 23 11 08    12 22  6 13; ...
     1938     3 21  6 40     6 22  1 59     9 23 16 56    12 22 12 03; ...
     1939     3 21 12 29     6 22  7 48     9 23 22 45    12 22 17 53; ...
     1940     3 20 18 18     6 21 13 36     9 23  4 34    12 21 23 42; ...
     1941     3 21  0 07     6 21 19 24     9 23 10 22    12 22  5 32; ...
     1942     3 21  5 56     6 22  1 12     9 23 16 11    12 22 11 22; ...
     1943     3 21 11 46     6 22  7 00     9 23 22 00    12 22 17 12; ...
     1944     3 20 17 35     6 21 12 48     9 23  3 49    12 21 23 01; ...
     1945     3 20 23 24     6 21 18 36     9 23  9 37    12 22  4 51; ...
     1946     3 21  5 13     6 22  0 25     9 23 15 26    12 22 10 41; ...
     1947     3 21 11 02     6 22  6 13     9 23 21 15    12 22 16 31; ...
     1948     3 20 16 52     6 21 12 01     9 23  3 03    12 21 22 20; ...
     1949     3 20 22 41     6 21 17 49     9 23  8 52    12 22  4 10; ...
     1950     3 21  4 30     6 21 23 37     9 23 14 41    12 22 10 00; ...
     1951     3 21 10 19     6 22  5 25     9 23 20 30    12 22 15 50; ...
     1952     3 20 16 08     6 21 11 13     9 23  2 18    12 21 21 39; ...
     1953     3 20 21 58     6 21 17 01     9 23  8 07    12 22  3 29; ...
     1954     3 21  3 47     6 21 22 50     9 23 13 56    12 22  9 19; ...
     1955     3 21  9 36     6 22  4 38     9 23 19 44    12 22 15 09; ...
     1956     3 20 15 25     6 21 10 26     9 23  1 33    12 21 20 58; ...
     1957     3 20 21 14     6 21 16 14     9 23  7 22    12 22  2 48; ...
     1958     3 21  3 04     6 21 22 02     9 23 13 11    12 22  8 38; ...
     1959     3 21  8 53     6 22  3 50     9 23 18 59    12 22 14 27; ...
     1960     3 20 14 42     6 21  9 38     9 23  0 48    12 21 20 17; ...
     1961     3 20 20 31     6 21 15 27     9 23  6 37    12 22  2 07; ...
     1962     3 21  2 20     6 21 21 15     9 23 12 25    12 22  7 57; ...
     1963     3 21  8 10     6 22  3 03     9 23 18 14    12 22 13 46; ...
     1964     3 20 13 59     6 21  8 51     9 23  0 03    12 21 19 36; ...
     1965     3 20 19 48     6 21 14 39     9 23  5 52    12 22  1 26; ...
     1966     3 21  1 37     6 21 20 27     9 23 11 40    12 22  7 16; ...
     1967     3 21  7 26     6 22  2 15     9 23 17 29    12 22 13 05; ...
     1968     3 20 13 16     6 21  8 03     9 22 23 18    12 21 18 55; ...
     1969     3 20 19 05     6 21 13 52     9 23  5 06    12 22  0 45; ...
     1970     3 21  0 54     6 21 19 40     9 23 10 55    12 22  6 35; ...
     1971     3 21  6 43     6 22  1 28     9 23 16 44    12 22 12 24; ...
     1972     3 20 12 32     6 21  7 16     9 22 22 33    12 21 18 14; ...
     1973     3 20 18 22     6 21 13 04     9 23  4 21    12 22  0 04; ...
     1974     3 21  0 11     6 21 18 52     9 23 10 10    12 22  5 54; ...
     1975     3 21  6 00     6 22  0 40     9 23 15 59    12 22 11 43; ...
     1976     3 20 11 49     6 21  6 29     9 22 21 47    12 21 17 33; ...
     1977     3 20 17 38     6 21 12 17     9 23  3 36    12 21 23 23; ...
     1978     3 20 23 28     6 21 18 05     9 23  9 25    12 22  5 13; ...
     1979     3 21  5 17     6 21 23 53     9 23 15 14    12 22 11 02; ...
     1980     3 20 11 06     6 21  5 41     9 22 21 02    12 21 16 52; ...
     1981     3 20 16 55     6 21 11 29     9 23  2 51    12 21 22 42; ...
     1982     3 20 22 44     6 21 17 17     9 23  8 40    12 22  4 31; ...
     1983     3 21  4 34     6 21 23 05     9 23 14 28    12 22 10 21; ...
     1984     3 20 10 23     6 21  4 54     9 22 20 17    12 21 16 11; ...
     1985     3 20 16 12     6 21 10 42     9 23  2 06    12 21 22 01; ...
     1986     3 20 22 01     6 21 16 30     9 23  7 54    12 22  3 50; ...
     1987     3 21  3 50     6 21 22 18     9 23 13 43    12 22  9 40; ...
     1988     3 20  9 40     6 21  4 06     9 22 19 32    12 21 15 30; ...
     1989     3 20 15 29     6 21  9 54     9 23  1 21    12 21 21 20; ...
     1990     3 20 21 18     6 21 15 42     9 23  7 09    12 22  3 09; ...
     1991     3 21  3 07     6 21 21 31     9 23 12 58    12 22  8 59; ...
     1992     3 20  8 56     6 21  3 19     9 22 18 47    12 21 14 49; ...
     1993     3 20 14 46     6 21  9 07     9 23  0 35    12 21 20 39; ...
     1994     3 20 20 35     6 21 14 55     9 23  6 24    12 22  2 28; ...
     1995     3 21  2 24     6 21 20 43     9 23 12 13    12 22  8 18; ...
     1996     3 20  8 13     6 21  2 31     9 22 18 01    12 21 14 08; ...
     1997     3 20 14 02     6 21  8 19     9 22 23 50    12 21 19 57; ...
     1998     3 20 19 52     6 21 14 07     9 23  5 39    12 22  1 47; ...
     1999     3 21  1 41     6 21 19 56     9 23 11 28    12 22  7 37; ...
     2000     3 20  7 30     6 21  1 44     9 22 17 16    12 21 13 27; ...
     2001     3 20 13 19     6 21  7 32     9 22 23 05    12 21 19 16; ...
     2002     3 20 19 08     6 21 13 20     9 23  4 54    12 22  1 06; ...
     2003     3 21  0 58     6 21 19 08     9 23 10 42    12 22  6 56; ...
     2004     3 20  6 47     6 21  0 56     9 22 16 31    12 21 12 46; ...
     2005     3 20 12 36     6 21  6 44     9 22 22 20    12 21 18 35; ...
     2006     3 20 18 25     6 21 12 32     9 23  4 08    12 22  0 25; ...
     2007     3 21  0 14     6 21 18 21     9 23  9 57    12 22  6 15; ...
     2008     3 20  6 04     6 21  0 09     9 22 15 46    12 21 12 04; ...
     2009     3 20 11 53     6 21  5 57     9 22 21 34    12 21 17 54; ...
     2010     3 20 17 42     6 21 11 45     9 23  3 23    12 21 23 44; ...
     2011     3 20 23 31     6 21 17 33     9 23  9 12    12 22  5 34; ...
     2012     3 20  5 20     6 20 23 21     9 22 15 01    12 21 11 23; ...
     2013     3 20 11 10     6 21  5 09     9 22 20 49    12 21 17 13; ...
     2014     3 20 16 59     6 21 10 57     9 23  2 38    12 21 23 03; ...
     2015     3 20 22 48     6 21 16 46     9 23  8 27    12 22  4 53; ...
     2016     3 20  4 37     6 20 22 34     9 22 14 15    12 21 10 42; ...
     2017     3 20 10 26     6 21  4 22     9 22 20 04    12 21 16 32; ...
     2018     3 20 16 16     6 21 10 10     9 23  1 53    12 21 22 22; ...
     2019     3 20 22 05     6 21 15 58     9 23  7 41    12 22  4 11; ...
     2020     3 20  3 54     6 20 21 46     9 22 13 30    12 21 10 01; ...
     2021     3 20  9 43     6 21  3 34     9 22 19 19    12 21 15 51; ...
     2022     3 20 15 32     6 21  9 23     9 23  1 07    12 21 21 41; ...
     2023     3 20 21 22     6 21 15 11     9 23  6 56    12 22  3 30; ...
     2024     3 20  3 11     6 20 20 59     9 22 12 45    12 21  9 20; ...
     2025     3 20  9 00     6 21  2 47     9 22 18 33    12 21 15 10; ...
     2026     3 20 14 49     6 21  8 35     9 23  0 22    12 21 20 59; ...
     2027     3 20 20 38     6 21 14 23     9 23  6 11    12 22  2 49; ...
     2028     3 20  2 28     6 20 20 11     9 22 11 59    12 21  8 39; ...
     2029     3 20  8 17     6 21  1 59     9 22 17 48    12 21 14 29; ...
     2030     3 20 14 06     6 21  7 48     9 22 23 37    12 21 20 18; ...
     2031     3 20 19 55     6 21 13 36     9 23  5 26    12 22  2 08; ...
     2032     3 20  1 44     6 20 19 24     9 22 11 14    12 21  7 58; ...
     2033     3 20  7 34     6 21  1 12     9 22 17 03    12 21 13 48; ...
     2034     3 20 13 23     6 21  7 00     9 22 22 52    12 21 19 37; ...
     2035     3 20 19 12     6 21 12 48     9 23  4 40    12 22  1 27; ...
     2036     3 20  1 01     6 20 18 36     9 22 10 29    12 21  7 17; ...
     2037     3 20  6 50     6 21  0 24     9 22 16 18    12 21 13 06; ...
     2038     3 20 12 40     6 21  6 13     9 22 22 06    12 21 18 56; ...
     2039     3 20 18 29     6 21 12 01     9 23  3 55    12 22  0 46; ...
     2040     3 20  0 18     6 20 17 49     9 22  9 44    12 21  6 36; ...
     2041     3 20  6 07     6 20 23 37     9 22 15 32    12 21 12 25; ...
     2042     3 20 11 56     6 21  5 25     9 22 21 21    12 21 18 15; ...
     2043     3 20 17 46     6 21 11 13     9 23  3 10    12 22  0 05; ...
     2044     3 19 23 35     6 20 17 01     9 22  8 58    12 21  5 54; ...
     2045     3 20  5 24     6 20 22 49     9 22 14 47    12 21 11 44; ...
     2046     3 20 11 13     6 21  4 37     9 22 20 36    12 21 17 34; ...
     2047     3 20 17 02     6 21 10 26     9 23  2 24    12 21 23 24; ...
     2048     3 19 22 52     6 20 16 14     9 22  8 13    12 21  5 13; ...
     2049     3 20  4 41     6 20 22 02     9 22 14 02    12 21 11 03; ...
     2050     3 20 10 30     6 21  3 50     9 22 19 50    12 21 16 53; ...
     2051     3 20 16 19     6 21  9 38     9 23  1 39    12 21 22 42; ...
     2052     3 19 22 08     6 20 15 26     9 22  7 28    12 21  4 32; ...
     2053     3 20  3 58     6 20 21 14     9 22 13 16    12 21 10 22; ...
     2054     3 20  9 47     6 21  3 02     9 22 19 05    12 21 16 12; ...
     2055     3 20 15 36     6 21  8 51     9 23  0 54    12 21 22 01; ...
     2056     3 19 21 25     6 20 14 39     9 22  6 42    12 21  3 51; ...
     2057     3 20  3 14     6 20 20 27     9 22 12 31    12 21  9 41; ...
     2058     3 20  9 04     6 21  2 15     9 22 18 20    12 21 15 30; ...
     2059     3 20 14 53     6 21  8 03     9 23  0 08    12 21 21 20; ...
     2060     3 19 20 42     6 20 13 51     9 22  5 57    12 21  3 10; ...
     2061     3 20  2 31     6 20 19 39     9 22 11 46    12 21  9 00; ...
     2062     3 20  8 20     6 21  1 27     9 22 17 34    12 21 14 49; ...
     2063     3 20 14 10     6 21  7 16     9 22 23 23    12 21 20 39; ...
     2064     3 19 19 59     6 20 13 04     9 22  5 12    12 21  2 29; ...
     2065     3 20  1 48     6 20 18 52     9 22 11 00    12 21  8 18; ...
     2066     3 20  7 37     6 21  0 40     9 22 16 49    12 21 14 08; ...
     2067     3 20 13 26     6 21  6 28     9 22 22 38    12 21 19 58; ...
     2068     3 19 19 16     6 20 12 16     9 22  4 26    12 21  1 47; ...
     2069     3 20  1 05     6 20 18 04     9 22 10 15    12 21  7 37; ...
     2070     3 20  6 54     6 20 23 52     9 22 16 04    12 21 13 27; ...
     2071     3 20 12 43     6 21  5 41     9 22 21 52    12 21 19 17; ...
     2072     3 19 18 32     6 20 11 29     9 22  3 41    12 21  1 06; ...
     2073     3 20  0 22     6 20 17 17     9 22  9 30    12 21  6 56; ...
     2074     3 20  6 11     6 20 23 05     9 22 15 18    12 21 12 46; ...
     2075     3 20 12 00     6 21  4 53     9 22 21 07    12 21 18 35; ...
     2076     3 19 17 49     6 20 10 41     9 22  2 56    12 21  0 25; ...
     2077     3 19 23 38     6 20 16 29     9 22  8 44    12 21  6 15; ...
     2078     3 20  5 28     6 20 22 17     9 22 14 33    12 21 12 05; ...
     2079     3 20 11 17     6 21  4 05     9 22 20 22    12 21 17 54; ...
     2080     3 19 17 06     6 20  9 54     9 22  2 10    12 20 23 44; ...
     2081     3 19 22 55     6 20 15 42     9 22  7 59    12 21  5 34; ...
     2082     3 20  4 44     6 20 21 30     9 22 13 48    12 21 11 23; ...
     2083     3 20 10 34     6 21  3 18     9 22 19 36    12 21 17 13; ...
     2084     3 19 16 23     6 20  9 06     9 22  1 25    12 20 23 03; ...
     2085     3 19 22 12     6 20 14 54     9 22  7 14    12 21  4 52; ...
     2086     3 20  4 01     6 20 20 42     9 22 13 02    12 21 10 42; ...
     2087     3 20  9 50     6 21  2 30     9 22 18 51    12 21 16 32; ...
     2088     3 19 15 40     6 20  8 19     9 22  0 39    12 20 22 22; ...
     2089     3 19 21 29     6 20 14 07     9 22  6 28    12 21  4 11; ...
     2090     3 20  3 18     6 20 19 55     9 22 12 17    12 21 10 01; ...
     2091     3 20  9 07     6 21  1 43     9 22 18 05    12 21 15 51; ...
     2092     3 19 14 56     6 20  7 31     9 21 23 54    12 20 21 40; ...
     2093     3 19 20 46     6 20 13 19     9 22  5 43    12 21  3 30; ...
     2094     3 20  2 35     6 20 19 07     9 22 11 31    12 21  9 20; ...
     2095     3 20  8 24     6 21  0 55     9 22 17 20    12 21 15 09; ...
     2096     3 19 14 13     6 20  6 43     9 21 23 09    12 20 20 59; ...
     2097     3 19 20 02     6 20 12 32     9 22  4 57    12 21  2 49; ...
     2098     3 20  1 52     6 20 18 20     9 22 10 46    12 21  8 39; ...
     2099     3 20  7 41     6 21  0 08     9 22 16 35    12 21 14 28; ...
     2100     3 20 13 30     6 21  5 56     9 22 22 23    12 21 20 18;];
end