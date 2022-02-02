REFIT Smart Home dataset
---
AFFILIATIONS
Steven Firth, Tom Kane, Vanda Dimitriou, Tarek Hassan, Farid Fouchal, Michael Coleman and Lynda Webb
Building Energy Research Group, School of Civil and Building Engineering, Loughborough University, UK
---
INTRO
This dataset is made publically available to support research in the area of end-use energy demand in the domestic building sector. 
Possible topics which this dataset might inform are:
1. Internal environmental conditions and energy demands in UK homes.
2. Analytic techniques for Smart Meter and Smart Home data.
3. Occupant behaviour in UK homes.
4. Model inputs for building energy models.
Please download and make use of this dataset in your own research and we are happy to answer any questions about the dataset and its use.
---
CITATION INFORMATION
Firth, Steven; Kane, Tom; Dimitriou, Vanda; Hassan, Tarek; Fouchal, Farid; Coleman, Michael; Webb, Lynda (): REFIT Smart Home dataset. figshare. 
https://dx.doi.org/10.17028/rd.lboro.2070091
---
DESCRIPTION OF DATASET
The REFIT research project (www.refitsmarthomes.org) ran from May 2012 to October 2015 and was a collaboration between Loughborough University, the University of Strathclyde and the University of East Anglia. 
The project carried out a study of Smart Home technologies in 20 UK homes which included a series of sensor measurements, building surveys and household interviews.
There are four files available for download:
1.ReadMe.txt: A ReadMe file with more details of the study and the dataset.
2.RefitXMLSchema.xsd: The RefitXML schema file which describes the structure and content of the REFIT_BUILDING_SURVEY.xml file.
3.REFIT_BUILDING_SURVEY.xml: A RefitXML file containing information about the buildings and the sensors placed in the buildings. 
4.REFIT_TIME_SERIES_VALUES.csv: A csv file containing the measurements recorded by the sensors placed in the buildings.
---
SUMMARY OF DATA COLLECTION
The dataset is for 20 homes located near to the town of Loughborough in the East Midlands region of the UK.
A building survey was carried out at each home, collecting data on building geometry, construction materials, occuapancy and energy services.
Each home has a selection of the following sensors and devices installed:
- CurrentCost mains clamps, to measure household mains electrical power load (data available from Strathclyde University, see below).
- Replacement gas meters, to measure household mains gas consumption.
- Hobo pendant or Hobo U12 sensors to measure room air temperature, relative humidity and light level.
- iButton temperature sensors to measure radiator surface temperature.
- CurrentCost individual appliance monitors, to measure plug electrical power loads (data available from Strathclyde University, see below).
- RWE Smart Home devices including programmable thermostatic radiator valves, interior and exterior motion detectors, door and window opening sensors and smoke alarms.
- British Gas Hive programmable thermostats.
In addition climate data was collected at the Loughborough University campus weather station.
---
TIMELINE
September 2013 to February 2014: Building surveys were carried out and monitoring sensors were placed in the buildings at or shortly after this time.
June 2014 and October 2014: Smart Home devices were installed in the buildings.
April 2015: Data collection finished.
---
DATA STATISTICS
Number of homes: 20
Number of spaces (rooms): 389
Number of radiators: 252
Number of showers: 34
Number of appliances: 618
Number of light bulbs: 672
Number of fixed heaters: 19
Number of surfaces: 2237
Number of openings: 970
Number of sensors: 1,567
Number of variables recorded by sensors and devices: 2,457
Number of time series readings: 25,312,397
---
USING THE DATA FILES
'REFIT_BUILDING_SURVEY.xml'
- This file can be opened using Notepad, an internet browser or an XML editor such as Microsoft XML Notepad. 
- Data can be extracted from the xml file by opening it in Excel or by using programming code (such as Python using the lxml package or Visual Basic for Applications using the DOMDOCUMENT60 object in the MSXML6.0 reference library).
- The descriptions of the elements (i.e 'Building', 'Space') and attributes (i.e. 'id', 'startDateTime') in this file can be found in the 'RefitXMLSchema.xsd' schema, under the 'annotation' nodes associated with each element or attribute. 
'REFIT_TIME_SERIES_VALUES.csv'
- Each row in this file represents a 'TimeSeriesValue' element in the RefitXML schema.
- This data was not placed in the 'REFIT_BUILDING_SURVEY.xml' file as the xml file would become too large and unusable. Instead the data was placed in this comma separated variable file. 
- The data is linked to the 'REFIT_BUILDING_SURVEY.xml' file using the 'TimeSeriesVariable/@id' column, which gives the 'id' attribute value of the parent TimeSeriesVariable element.
- This file contains 25,312,397 rows and is too large to be opened by Notepad or Excel.
- It can be opened by creating a 'link table' in Microsoft Access or using another database programme, or the data can be extracted using programming code (such as Python using the pandas package).
---
ADDITIONAL REFIT DATASETS - ELECTRICITY DATA, UNIVERSITY OF STRATHCLYDE
Murray, David and Stankovic, Lina (2016). REFIT: Electrical Load Measurements (Cleaned) - available at: http://dx.doi.org/10.15129/9ab14b0e-19ac-4279-938f-27f643078cec
---
ADDITIONAL REFIT DATASETS - SURVEY DATA AND SMART HOME INTERVIEWS, UNIVERSITY OF EAST ANGLIA
Wilson, Charlie and Hargreaves, Tom (2016). REFIT: Personalised retrofit decision support tools for UK homes using smart home technology. Phase 1: Survey data. [Data Collection]. Colchester, Essex: UK Data Archive.  10.5255/UKDA-SN-852366 - available at: https://dx.doi.org/10.5255/UKDA-SN-852366
Wilson, Charlie and Hargreaves, Tom and Hauxwell-Baldwin, Richard (2016). REFIT: Personalised retrofit decision support tools for UK homes using smart home technology. Phase 2: Smart home interviews. [Data Collection]. Colchester, Essex: UK Data Archive.  10.5255/UKDA-SN-852367  - available at: https://dx.doi.org/10.5255/UKDA-SN-852367
---
FURTHER INFORMATION
For more information about the REFIT project please contact Dr Steven Firth: s.k.firth@lboro.ac.uk, 01509 228546, http://www.lboro.ac.uk/departments/civil-building/staff/firthsteven/
---
ACKNOWLEDGEMENTS
This work has been carried out as part of the REFIT project (‘Personalised Retrofit Decision Support Tools for UK Homes using Smart Home Technology’, Grant Reference EP/K002457/1). 
REFIT is a consortium of three universities - Loughborough, Strathclyde and East Anglia - and ten industry stakeholders funded by the Engineering and Physical Sciences Research Council (EPSRC) under the Transforming Energy Demand in Buildings through Digital Innovation (BuildTEDDI) funding programme. 
For more information see: www.epsrc.ac.uk and www.refitsmarthomes.org

