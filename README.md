![alt text](https://nhsbtdbe.blob.core.windows.net/umbraco-assets-corp/25124/blood-g4e28dcb97_1920.jpg)


# WBA_Semi-automatized_EC50-calculation
Human Whole blood assay (WBA) -  blood was treated with different compounds using different concentrations in order to determine EC50 (Half maximal effective concentration) values regarding different cytokines.


  1. Human blood was equally distributed onto a 96-well plate (225µl/well)
  
  2. Prediluted compounds (A and B in triplicates) were added to the wells respectively (25µl/well)
  
  3. After 24h incubation, plate was centrifuged and the serum part was collected (~50µl/well)
  
  4. Collected supernatant was used for cytokine measurment in an MSD system (in each well 10 different spots are ready to measure individual analyte concentration (IFNb, IL-10...) symulteniously)
  
  5. The titration allows to determine EC50 values regarding all secreted cytokines




## Test
The calculation was performed in GraphPad Prism in order to double check the accuracy 

### Needed files:
1. EC50_test.R
2. EC50_training_data.xlsx
3. Test_calculation_GraphPad_2_.JPG
4. Test_calculation_GraphPad_2_.JPG

## Experiment
Experiment data was used to compare two different commercially available STING antagonist ([A](https://www.tocris.com/products/sting-agonist-c53_7741),  and [B](https://www.selleckchem.com/products/diabzi-sting-agonist-compound-3.html)) to determine capability to induce cytokine secretion then calculate EC50 values for each cytokine respecively. 

### Layout: triplicates (concentration is in µM)

![alt text](https://github.com/AdamAdonyi/WBA_Semiautomatized_EC50-calculation/blob/main/Layout.JPG)

### Investigated cytokines: IL-10, IL-6, IFNg, IFNa, IP-10, MCP-1, IL-2, GRZMb, IFNb, TNFa

### Needed files:
1. 2BM2HAMM84
2. WBA.R
3. Result_EC50_A.JPG
4. Result_EC50_B.JPG

# Results
![alt text](https://github.com/AdamAdonyi/WBA_Semiautomatized_EC50-calculation/blob/main/Result_EC50_A.JPG) |
![alt text](https://github.com/AdamAdonyi/WBA_Semiautomatized_EC50-calculation/blob/main/Result_EC50_B.JPG)

# References:

[Used MSD technology](https://www.mesoscale.com/en/products_and_services/assay_kits/multiplex_assay_kits/)
