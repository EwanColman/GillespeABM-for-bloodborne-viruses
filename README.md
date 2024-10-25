This is the Repo for my implementation of an agent-based model that uses the Gillespe algorithm. It's for a project on transmissible bloodborne diseases in people who inject drugs.
![HIV](https://github.com/user-attachments/assets/7921586d-8df0-4740-825a-c9337615d5b0)
![HCV](https://github.com/user-attachments/assets/3d3391d7-9e9c-47a9-aa21-c5f8fef70da3)

Read the code explainer written in Jupyter notebook

Other scripts are in development for a specific project. 
HIV_and_HCV_model_with_pop_dynamics_and_interventions.py performs one run of the model. HIV_and_HCV_model_with_pop_dynamics_and_interventions_many_runs.py runs the model lots of times and saves the output as a pickle file. 
plot_intervention_comparison.py requires pickle files (created by the script above) to exist. Various parts of the code should be commented or uncommented to produce the three necessary outputs (corresponding to different intervention strategies).
