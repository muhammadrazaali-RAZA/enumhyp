package at.univie.federated.server;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

/**
 * Handles file operations for saving received client data.
 * Simple, focused responsibility: save data to files.
 */
public class FileManager {
    
    private static final String OUTPUT_DIR = "serverReceived";
    private final ErrorHandler errorHandler;
    
    public FileManager(ErrorHandler errorHandler) {
        this.errorHandler = errorHandler;
    }
    
    /**
     * Saves data lines to a file named server{nodeId}.txt in the output directory.
     * 
     * @param nodeId The node identifier
     * @param dataLines List of data lines to save
     * @return Path to the saved file, or null if failed
     */
    public Path saveData(String nodeId, List<String> dataLines) {
        if (nodeId == null || dataLines == null || dataLines.isEmpty()) {
            errorHandler.logWarning("Cannot save data: nodeId or dataLines is null/empty");
            return null;
        }
        
        try {
            Path dir = Paths.get(OUTPUT_DIR);
            if (!Files.exists(dir)) {
                Files.createDirectories(dir);
            }
            
            Path outputFile = dir.resolve("server" + nodeId + ".txt");
            Files.write(outputFile, dataLines);
            
            return outputFile;
        } catch (IOException e) {
            errorHandler.handleFileError(e, "server" + nodeId + ".txt", "saveData");
            return null;
        }
    }
}

