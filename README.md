# Behav3d Docker Installation Guide

This guide will walk you through the steps to install Docker and run Behav3d using a Docker container.

## Installing Docker

To install Docker, follow the instructions provided on the [Docker website](https://www.docker.com/get-started). Docker provides comprehensive installation guides for various operating systems.

## Running Behav3d Docker Container

Once Docker is installed, make sure that Docker is open and then you can run Behav3d using the following command:

```bash
docker run -p 8888:8888 -v [C:/Users/Malieva/Documents/Roca/docker]:/app/data imaigenelab/behav3d:1.0
```

Important: Replace `[C:/Users/Malieva/Documents/Roca/docker]` with the path where your data is stored. Make sure to provide the correct path to access your data and then remove the brackets.

ALL the data used in the pipeline must be in the parent folder or in a subfolder

After launching the command, you can access Behav3d by opening any web browser and typing the following URL in the address bar:

```
localhost:8888
```

This will open the Behav3d application in your browser, allowing you to start working with your data.

## Stopping and Restarting Behav3d Docker Container
To stop the Behav3d Docker container, follow these steps:

Open Docker Desktop.
1. Navigate to the Containers section. You can find it by clicking on the Docker icon in the system tray and selecting "Containers".
2. In the Containers section, you will see a running containers. Locate the Behav3d container.
3. Click on the "Stop" button next to the Behav3d container to stop it.
4. To restart the Behav3d Docker container, simply click on the "Start" button next to the stopped container.
## Additional Notes

- ## FIRST TIME USING DO NOT CLOSE TERMINAL THAT LAUNCHED THE COMMAND
- Ensure that Docker is running before executing the command.
- Make sure to replace the placeholder path with the actual path to your data directory.
- You may need appropriate permissions to access certain directories.
- For more information on Behav3d and its usage, refer to the official documentation.
