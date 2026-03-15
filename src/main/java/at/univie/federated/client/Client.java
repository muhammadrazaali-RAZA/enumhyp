package at.univie.federated.client;

import java.io.*;
import java.net.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;


/**
 * Client class: Entry point for the federated client application.
 *
 * This client establishes a TCP connection to a server running on localhost:1234.
 * It reads the UCC results from abalone_p1.txt and sends the data to the server.
 *
 * Key Features:
 * - Socket-based communication using java.net.Socket
 * - Reads and transmits file data (abalone_p1.txt)
 * - Buffered I/O streams for efficient data transfer
 * - Graceful resource cleanup in finally block
 */
public class Client {

    public void activateClient(String fileName, String nodeId) {
        System.out.println("Socket Client Started for node: " + nodeId + " with file: " + fileName);

        String host = System.getenv().getOrDefault("SERVER_HOST", "localhost");
        int port = Integer.parseInt(System.getenv().getOrDefault("SERVER_PORT", "1234"));
        
        // Get file path from environment or use default
        String dataDir = System.getenv().getOrDefault("DATA_DIR", "data");
        Path filePath = Paths.get(dataDir, fileName);

        Socket socket = null;
        InputStreamReader inputStreamReader = null;
        OutputStreamWriter outputStreamWriter = null;
        BufferedReader bufferedReader = null;
        BufferedWriter bufferedWriter = null;

        try {
            socket = new Socket(host, port);
            inputStreamReader = new InputStreamReader(socket.getInputStream());
            outputStreamWriter = new OutputStreamWriter(socket.getOutputStream());
            bufferedReader = new BufferedReader(inputStreamReader);
            bufferedWriter = new BufferedWriter(outputStreamWriter);

            // 1) Send a greeting with node ID
            String greeting = "HELLO from " + nodeId;
            bufferedWriter.write(greeting);
            bufferedWriter.newLine();
            bufferedWriter.flush();
            System.out.println("Sent: " + greeting);
            String serverResponse = bufferedReader.readLine();
            System.out.println("Server: " + serverResponse);

            // 2) Read file and send data
            if (Files.exists(filePath)) {
                System.out.println("Reading file: " + filePath);
                List<String> fileLines = Files.readAllLines(filePath);
                
                // Send file name and line count
                bufferedWriter.write("FILE_DATA");
                bufferedWriter.newLine();
                bufferedWriter.write(fileName);
                bufferedWriter.newLine();
                bufferedWriter.write(String.valueOf(fileLines.size()));
                bufferedWriter.newLine();
                bufferedWriter.flush();
                
                // Send file content line by line with nodeId attached
                for (String line : fileLines) {
                    // Attach nodeId to each line in format: nodeId|data
                    bufferedWriter.write(nodeId + "|" + line);
                    bufferedWriter.newLine();
                }
                bufferedWriter.flush();
                
                System.out.println("Sent file data: " + fileLines.size() + " lines");
                // Wait for server acknowledgment
                serverResponse = bufferedReader.readLine();
                System.out.println("Server: " + serverResponse);
            } else {
                System.err.println("File not found: " + filePath);
                bufferedWriter.write("FILE_NOT_FOUND");
                bufferedWriter.newLine();
                bufferedWriter.flush();
            }

            // 3) Send QUIT to close the session cleanly
            String quitMsg = "QUIT";
            bufferedWriter.write(quitMsg);
            bufferedWriter.newLine();
            bufferedWriter.flush();
            System.out.println("Sent: " + quitMsg);
            serverResponse = bufferedReader.readLine();
            System.out.println("Server: " + serverResponse);

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (socket != null) socket.close();
                if (inputStreamReader != null) inputStreamReader.close();
                if (outputStreamWriter != null) outputStreamWriter.close();
                if (bufferedReader != null) bufferedReader.close();
                if (bufferedWriter != null) bufferedWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}