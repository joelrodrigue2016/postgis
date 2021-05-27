select *
from public.nyc_yellow_taxi_trips_2016_06_01
limit 1000;


-- Add column
ALTER TABLE nyc_yellow_taxi_trips_2016_06_01 ADD COLUMN geom_pickup geography(POINT,4326);

-- Now fill that column with the lat/long
UPDATE nyc_yellow_taxi_trips_2016_06_01
SET geom_pickup = ST_SetSRID(
                            ST_MakePoint(pickup_longitude,pickup_latitude),4326
                           )::geography;

-- Add a GiST index
CREATE INDEX yc_yellow_taxi_trips_2016_06_01_pts_idx ON yc_yellow_taxi_trips_2016_06_01 USING GIST (geom_pickup);
