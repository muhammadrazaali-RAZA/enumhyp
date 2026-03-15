package at.univie.federated.server;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Handles different types of messages from clients.
 * Simple, focused responsibility: process and respond to messages.
 */
public class MessageHandler {
    
    private final FileManager fileManager;
    private final ErrorHandler errorHandler;
    
    public MessageHandler(FileManager fileManager, ErrorHandler errorHandler) {
        this.fileManager = fileManager;
        this.errorHandler = errorHandler;
    }
    
    /**
     * Handles a HELLO message from client.
     */
    public void handleHello(String message, BufferedWriter writer) throws IOException {
        System.out.println("Client: " + message);
        sendResponse(writer, "Hello received. Ready for data.");
    }
    
    /**
     * Handles FILE_DATA message - receives and saves file data.
     */
    public void handleFileData(BufferedReader reader, BufferedWriter writer) throws IOException {
        try {
            String fileName = reader.readLine();
            if (fileName == null) {
                errorHandler.logWarning("Received null filename");
                return;
            }
            
            String lineCountStr = reader.readLine();
            if (lineCountStr == null) {
                errorHandler.logWarning("Received null line count");
                return;
            }
            
            int lineCount;
            try {
                lineCount = Integer.parseInt(lineCountStr);
            } catch (NumberFormatException e) {
                errorHandler.handleNumberFormatException(e, "handleFileData", lineCountStr);
                sendResponse(writer, "Error: Invalid line count format");
                return;
            }
            
            System.out.println("Receiving file: " + fileName);
            System.out.println("Expected lines: " + lineCount);
            
            String nodeId = null;
            List<String> dataLines = new ArrayList<>();
            
            // Read all lines and extract nodeId from first line
            for (int i = 0; i < lineCount; i++) {
                String line = reader.readLine();
                if (line != null) {
                    if (line.contains("|")) {
                        String[] parts = line.split("\\|", 2);
                        // Extract nodeId from first line
                        if (nodeId == null) {
                            nodeId = parts[0];
                        }
                        String data = parts.length > 1 ? parts[1] : "";
                        dataLines.add(data);
                    } else {
                        dataLines.add(line);
                    }
                }
            }
            
            printReceivedData(fileName, dataLines, nodeId);
            saveDataToFile(nodeId, dataLines);
            
            sendResponse(writer, "File data received successfully. Total lines: " + dataLines.size());
        } catch (IOException e) {
            errorHandler.handleIOException(e, "handleFileData");
            throw e; // Re-throw to let caller handle
        }
    }
    
    /**
     * Handles FILE_NOT_FOUND message.
     */
    public void handleFileNotFound(BufferedWriter writer) throws IOException {
        System.out.println("Client: File not found on client side.");
        sendResponse(writer, "File not found acknowledged.");
    }
    
    /**
     * Handles QUIT message.
     */
    public boolean handleQuit(BufferedWriter writer) throws IOException {
        System.out.println("Client requested quit.");
        sendResponse(writer, "Session terminated. Goodbye.");
        return true; // Signal to break the loop
    }
    
    /**
     * Handles unknown messages.
     */
    public void handleUnknown(String message, BufferedWriter writer) throws IOException {
        System.out.println("Client: " + message);
        sendResponse(writer, "Message Received.");
    }
    
    // Private helper methods
    
    private void printReceivedData(String fileName, List<String> dataLines, String nodeId) {
        System.out.println("Received file data from client:");
        System.out.println("--- File: " + fileName + " ---");
        for (String line : dataLines) {
            if (nodeId != null) {
                System.out.println("[" + nodeId + "] " + line);
            } else {
                System.out.println(line);
            }
        }
        System.out.println("--- End of " + fileName + " ---");
    }
    
    private void saveDataToFile(String nodeId, List<String> dataLines) {
        if (nodeId != null && !dataLines.isEmpty()) {
            java.nio.file.Path savedFile = fileManager.saveData(nodeId, dataLines);
            if (savedFile != null) {
                System.out.println("Data saved to: " + savedFile);
            }
        }
    }
    
    private void sendResponse(BufferedWriter writer, String message) throws IOException {
        try {
            writer.write(message);
            writer.newLine();
            writer.flush();
        } catch (IOException e) {
            errorHandler.handleIOException(e, "sendResponse");
            throw e;
        }
    }
}
