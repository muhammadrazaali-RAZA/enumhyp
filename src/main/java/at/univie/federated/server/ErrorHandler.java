package at.univie.federated.server;

import java.io.IOException;
import java.net.SocketException;
import java.net.SocketTimeoutException;

/**
 * Centralized error handling for the server.
 * Simple responsibility: handle and log all exceptions consistently.
 */
public class ErrorHandler {
    
    /**
     * Handles general exceptions with appropriate logging.
     */
    public void handleException(Exception e, String context) {
        String errorMessage = "Error in " + context + ": " + e.getMessage();
        System.err.println(errorMessage);
        e.printStackTrace();
    }
    
    /**
     * Handles IOException - typically file or network I/O errors.
     */
    public void handleIOException(IOException e, String context) {
        String errorMessage = "I/O Error in " + context + ": " + e.getMessage();
        System.err.println(errorMessage);
        e.printStackTrace();
    }
    
    /**
     * Handles socket-related exceptions.
     */
    public void handleSocketException(SocketException e, String context) {
        String errorMessage = "Socket Error in " + context + ": " + e.getMessage();
        System.err.println(errorMessage);
        System.err.println("Connection may have been closed unexpectedly.");
        e.printStackTrace();
    }
    
    /**
     * Handles socket timeout exceptions.
     */
    public void handleSocketTimeout(SocketTimeoutException e, String context) {
        String errorMessage = "Socket Timeout in " + context + ": " + e.getMessage();
        System.err.println(errorMessage);
        System.err.println("Client connection timed out.");
        e.printStackTrace();
    }
    
    /**
     * Handles number format exceptions (e.g., parsing line count).
     */
    public void handleNumberFormatException(NumberFormatException e, String context, String invalidValue) {
        String errorMessage = "Number Format Error in " + context + ": Cannot parse '" + invalidValue + "'";
        System.err.println(errorMessage);
        System.err.println("Expected a valid number but received: " + invalidValue);
        e.printStackTrace();
    }
    
    /**
     * Handles file operation errors.
     */
    public void handleFileError(IOException e, String fileName, String operation) {
        String errorMessage = "File Error during " + operation + " for file: " + fileName;
        System.err.println(errorMessage);
        System.err.println("Error: " + e.getMessage());
        e.printStackTrace();
    }
    
    /**
     * Handles connection errors.
     */
    public void handleConnectionError(Exception e, String clientInfo) {
        String errorMessage = "Connection Error with client: " + clientInfo;
        System.err.println(errorMessage);
        System.err.println("Error: " + e.getMessage());
        e.printStackTrace();
    }
    
    /**
     * Handles resource cleanup errors.
     */
    public void handleCleanupError(IOException e, String resource) {
        String errorMessage = "Error closing resource: " + resource;
        System.err.println(errorMessage);
        System.err.println("Warning: Resource may not have been properly closed.");
        e.printStackTrace();
    }
    
    /**
     * Handles null pointer exceptions with context.
     */
    public void handleNullPointer(NullPointerException e, String context) {
        String errorMessage = "Null Pointer Error in " + context + ": A required object was null";
        System.err.println(errorMessage);
        System.err.println("Please check that all required objects are properly initialized.");
        e.printStackTrace();
    }
    
    /**
     * Logs a warning message.
     */
    public void logWarning(String message) {
        System.err.println("WARNING: " + message);
    }
    
    /**
     * Logs an info message about errors (non-critical).
     */
    public void logErrorInfo(String message) {
        System.err.println("INFO: " + message);
    }
}









