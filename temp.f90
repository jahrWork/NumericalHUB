   \newcommand{\twographs}[4]
     {
      \begin{figure}[htpb]
         \begin{minipage}[t]{0.5\textwidth} {#1} \end{minipage}
         \begin{minipage}[t]{0.5\textwidth} {#2} \end{minipage}
         \caption{#3} \label{#4} 
      \end{figure}  
     } 
                           
       \twographs
        {\input{./results/myexampleDa.tex}}
        {\input{./results/myexampleDa.tex}}
        {Heinon-Heiles system solution. 
        (a) Trajectory of the star $(x,y)$. 
        (b) Projection $(x,\dot{x})$ of the solution. 
        }
        {fig:exampleDa} 
