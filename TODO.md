# TODO:

- Wait until paper_data finished running and combine the following files using debug.ipynb:
    - data/20220314/paper_data_20220314_0.csv
    - data/20220314/paper_data_20220314_1.csv
    - data/20220314/paper_data_20220314_3.csv
    - data/paper_data_20220314.csv
    
- Expand definition of projects and groups.
- Re-create research group networks.


## OTHER

- Get departments from IGTP by using the information on their website.
- Check if isolated nodes are important. If they are we might need to recalculate edges based not only on joint publications but also affiliation to same department.


## NOTES:
Problem:
- Running the scraper from JupyterLab or directly from a script produces different results. 

Solution:
- This is actually because JupyterLab is inside in a Docker container. If you run the script from inside the docker container you get identical results.

