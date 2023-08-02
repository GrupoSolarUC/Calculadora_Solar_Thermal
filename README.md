# Calculadora_Solar_Thermal

Librería hecha para simular un sistema térmico solar como el que se ilustra en la imagen "Esquema_Sistema.jpg", en la carpeta "Files". El código está optimizado para obtener datos del "Explorador Solar" del Ministerio de Energía de Chile. Por lo tanto, está pensado para realizar estimaciones dentro del territorio nacional chileno. Para realizar una simulación es necesario seguir los siguientes pasos:
   
1. Instalar python y las demás librerías que se requieren. Las versiones recomendadas se especifican en la imagen "Librerias.jpg", en la carpeta "Files".

2. Descargar un archivo de datos climáticos (actualmente existe un archivo de ejemplo para evitar este paso; su nombre es DHTMY_E_7IH2CL.csv) siguiendo los pasos especificados en la imagen "Descargar_Archivos_Climaticos.jpg" que se encuentra en la carpeta "Files".

3. Llenar el archivo "Introducir_parametros.xlsx" con los parámetros del sistema y de la simulación (actualmente el archivo tiene parámetros de ejemplo).

4. Correr el Script "Simulate_System.py"

5. Esperar que termine de correr el código. Debiese haber una carpeta nueva con el nombre "Resultados_Simulacion_XX", donde XX corresponde al número de la simulación especificado por el usuario.

