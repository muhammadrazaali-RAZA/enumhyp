package at.univie.federated.server;

import java.io.*;
import java.net.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Server class: Entry point for the federated server application.
 * 
 * Simple responsibility: Accept connections and delegate to handlers.
 * Follows KISS principle - keeps it simple and focused.
 */
public class Server {

    private static final int PORT = 1234;
    private static final int THREAD_POOL_SIZE = 10;

    public void activateServer() throws IOException {
        System.out.println("Socket Server Started !");

        ErrorHandler errorHandler = new ErrorHandler();
        ServerSocket serverSocket = new ServerSocket(PORT);
        ExecutorService executorService = Executors.newFixedThreadPool(THREAD_POOL_SIZE);
        FileManager fileManager = new FileManager(errorHandler);
        MessageHandler messageHandler = new MessageHandler(fileManager, errorHandler);

        while (true) {
            try {
                Socket socket = serverSocket.accept();
                System.out.println("Client connected from: " + socket.getRemoteSocketAddress());
                
                executorService.submit(new ClientHandler(socket, messageHandler, errorHandler));
            } catch (IOException e) {
                errorHandler.handleIOException(e, "activateServer - accept connection");
            } catch (Exception e) {
                errorHandler.handleException(e, "activateServer");
            }
        }
    }
    
    /**
     * Handles individual client connections in separate threads.
     * Simple responsibility: Read messages and delegate to MessageHandler.
     */
    private static class ClientHandler implements Runnable {
        private final Socket socket;
        private final MessageHandler messageHandler;
        private final ErrorHandler errorHandler;
        
        public ClientHandler(Socket socket, MessageHandler messageHandler, ErrorHandler errorHandler) {
            this.socket = socket;
            this.messageHandler = messageHandler;
            this.errorHandler = errorHandler;
        }
        
        @Override
        public void run() {
            BufferedReader reader = null;
            BufferedWriter writer = null;

            try {
                reader = new BufferedReader(new InputStreamReader(socket.getInputStream()));
                writer = new BufferedWriter(new OutputStreamWriter(socket.getOutputStream()));

                while (true) {
                    String message = reader.readLine();

                    if (message == null) {
                        System.out.println("Client disconnected.");
                        break;
                    }

                    boolean shouldQuit = processMessage(message, reader, writer);
                    if (shouldQuit) {
                        break;
                    }
                }

            } catch (java.net.SocketException e) {
                errorHandler.handleSocketException(e, "ClientHandler.run");
            } catch (java.net.SocketTimeoutException e) {
                errorHandler.handleSocketTimeout(e, "ClientHandler.run");
            } catch (IOException e) {
                errorHandler.handleIOException(e, "ClientHandler.run");
            } catch (Exception e) {
                errorHandler.handleException(e, "ClientHandler.run");
            } finally {
                closeResources(reader, writer);
            }
        }
        
        private boolean processMessage(String message, BufferedReader reader, BufferedWriter writer) {
            try {
                if (message.startsWith("HELLO")) {
                    messageHandler.handleHello(message, writer);
                    return false;
                } else if (message.equals("FILE_DATA")) {
                    messageHandler.handleFileData(reader, writer);
                    return false;
                } else if (message.equals("FILE_NOT_FOUND")) {
                    messageHandler.handleFileNotFound(writer);
                    return false;
                } else if (message.equalsIgnoreCase("QUIT")) {
                    return messageHandler.handleQuit(writer);
                } else {
                    messageHandler.handleUnknown(message, writer);
                    return false;
                }
            } catch (IOException e) {
                errorHandler.handleIOException(e, "processMessage");
                return true; // Break loop on I/O error
            }
        }
        
        private void closeResources(BufferedReader reader, BufferedWriter writer) {
            try {
                if (socket != null) socket.close();
                if (reader != null) reader.close();
                if (writer != null) writer.close();
                System.out.println("Connection closed.");
            } catch (IOException e) {
                errorHandler.handleCleanupError(e, "socket/streams");
            }
        }
    }
}
