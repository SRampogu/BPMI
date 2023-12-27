import streamlit as st
import pandas as pd
#from PIL import Image
import subprocess
import os
import base64
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
import time

#-------------------------------------------------




#------------------------------------------------------------------------------------------
#st.title(':orange[Bioactivity Prediction of Plausible Mpro Inhibitors  (BPMI)]')
st.markdown("<h1 style='text-align: center; color: purple;'>Bioactivity Prediction of Plausible Mpro Inhibitors  (BPMI)</h1>", unsafe_allow_html=True)

def lipinski_check(SMILES):
    
    df=Chem.MolFromSmiles(SMILES)
    b=Descriptors.ExactMolWt(df)
    c=Descriptors.NumHAcceptors(df)
    d=Descriptors.NumHDonors(df)
    e=Descriptors.MolLogP(df)
    #df3=(b,c,d,e)
    df3=pd.Series({"MW":b, "nHA":c, "nHD":d,"LogP":e})
    #df4=pd.Series.to_frame(df3).T #(.T= transpose)
    return(df3)
#----------------------------------------------------------------------


#------------------------------------------------------------------------------------------------
# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')
#----------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
st.sidebar.header(":blue[Upload Data to Predict]:jigsaw:")
#uploaded_files = st.sidebar.file_uploader("Choose a CSV file", accept_multiple_files=False)
uploaded_file = st.sidebar.file_uploader("**:orange[Choose a file in CSV format]**", accept_multiple_files=False)
if uploaded_file is not None:
     #if st.sidebar.button('SUBMIT'):
         df0 = pd.read_csv(uploaded_file)
         #st.header("Data to Predict")
         st.header(' :orange[Data to Predict] ')
         st.write(df0)
         df0.to_csv('molecule.smi', sep = '\t', header = False, index = False)
         with st.spinner("In PROCESS!!..."):
                   desc_calc()
                   st.success("DONE!")       
          #Read in calculated descriptors and display the dataframe
         #st.header('**Calculated Molecular Descriptors**')
         st.header(' :orange[Calculated Molecular Descriptors] ')
         #st.write("**This data is used to predict the activity of compounds**")
         st.markdown("**:violet[*This data is used to predict the activity of compounds*]**")
         
         desc = pd.read_csv('descriptors_output.csv')
         st.write(desc)
         st.write(desc.shape)
       #---------------------------
        
        
        #---------------
         
         st.download_button("**:green[Download Calculated Descriptors]**:inbox_tray:",
                          desc.to_csv( ),
                          file_name = 'Descriptors.csv',
                          mime = 'text/csv')
#---------------------------------------------------------------------------------------
#Predicting the pIC50 values
         #st.header('Activity Prediction')
         st.header(' :orange[Activity Prediction] ')
         LipinModel = pickle.load(open('mixx_rfg.pkl', 'rb'))
         #st.write(LipinModel)
         df9=desc.drop(['Name'],axis=1)
         #st.write(df9)
         prediction = LipinModel.predict(df9)
         prediction_output = pd.Series(prediction, name='pIC50')
         df10=desc['Name']
         df11 = pd.concat([df10, prediction_output], axis=1)
         st.write('Predicted pIC50 values')
         st.write(df11)
         st.write(df11.shape)
         st.download_button('**:green[Download Predicted Activity]**:inbox_tray:',
                          df11.to_csv( ),
                          file_name = 'Predicted Activity.csv',
                          mime = 'text/csv')
#---------------------------------------------------------------------------------------------------
         #st.header ("Lipinski Prediction for a Dataframe")
         st.header(' :orange[Lipinski Prediction for a Dataframe] ')
         #df100=(df0.Name)
         

         
         #df100=lipinski_check(df)
         PandasTools.AddMoleculeColumnToFrame(df0, "SMILES")
         df200 = df0[~df0['ROMol'].isnull()]
         #st.write(df200)
         df150=df200.drop(['Name'], axis=1)
         #df200.to_csv('new_file.csv')
         df150.to_csv('new_file.csv')
         #st.write(df150)
         df400=df150["SMILES"].apply(lipinski_check)

         df20=df0.Name
         #df600=pd.concat([df20,df300], axis=1)
         df600=pd.concat([df20,df400], axis=1)
         
         st.write(df600)
         
         st.write(df600.shape)
         #st.write(df100)
        
         
         st.download_button('**:green[Download Predicted Lipinski]**:inbox_tray:',
                        
                          df600.to_csv( ),
                          file_name = 'Predicted Lipinski.csv',
                          mime = 'text/csv',
                          )
#--------------------------------------------------------------------------------------------
         
         st.sidebar.header('**:orange[2D structure visualisation]**')
         compound_smiles = st.sidebar.text_input('SMILES Please','COC1=C(C=CC(=C1)CC=C)O')
         m = Chem.MolFromSmiles(compound_smiles)
         im=Draw.MolToImage(m)
         st.sidebar.image(im)

#3D structure download
         m5 = Chem.MolFromSmiles(compound_smiles)
         m7 = Chem.AddHs(m5)
         m6= (Chem.MolToMolBlock(m7)) 
         st.sidebar.download_button(  label='3D.sdf', data=m6,file_name='3D.sdf') 

 #------------------------------------------------------------------------------------------------------        
#Structure download

         from io import BytesIO
         buf = BytesIO()
         im.save(buf, format="JPEG")
         byte_im = buf.getvalue()


         btn = st.sidebar.download_button(
                 label="2D Structure Image.png",
                 data=byte_im,
                 file_name="2D Structure Image.png",
                 mime="image/jpeg",
                  )
 
#-----------------------------------------------------------------------
#-------------------------------------------------
   
        

else:
      st.info('**:orange[Please Upload Data to Predict]**')
#-------------------------------------------------------   
