//package at.univie.federated;
//
//import at.univie.federated.aggregator.UccAggregator;
//import at.univie.federated.compared.FileComparator;
//import at.univie.federated.client.Client;
//import at.univie.enumhyp.GenerateEnumerate;
//import at.univie.federated.server.Server;
//
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.List;
//
//
//public class Main {
//    public static void main(String[] args) throws Exception {
//
//        boolean runServerClient = true;
//        boolean runEnumhyp = true;
//
//        // Variables for client management
//        // String data = "abalone";
//        String data = "fdReduced";
////        String data = "flight";
//
//        int length =   10;
//
//        boolean genGraph   = true;
//        boolean genEnuTxt  = true;
//
//        if (runEnumhyp) {
//                GenerateEnumerate runGenerateEnumerate = new GenerateEnumerate();
//                runGenerateEnumerate.generateEnumerate(data , length, genGraph, genEnuTxt);
//        }
//
//        if (runServerClient) {
//            // Start server in a separate thread
//
//            Server server = new Server();
//            Thread serverThread = new Thread(() -> {
//                try {
//                    server.activateServer();
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            });
//            serverThread.setDaemon(true);
//            serverThread.start();
//
//            // Wait a bit for server to start
//            Thread.sleep(1000);
//
//            // Start clients in separate threads to allow concurrent execution
//            List<Thread> clientThreads = new ArrayList<>();
//
//            for (int index = 0; index < length; index++) {
//                Client client = new Client();
//                final int currentPart = index + 1;
//                Thread clientThread = new Thread(() -> {
//                    client.activateClient(data + "/"+data + "_p" + currentPart + ".txt", "part_" + currentPart);
//                });
//                clientThreads.add(clientThread);
//                clientThread.start();
//            }
//
//            // Wait for all clients to finish
//            for (Thread thread : clientThreads) {
//                thread.join();
//            }
//
//            System.out.println("All clients completed. Server continues running...");
//        }
//
//        // Run aggregation of the UCCs
////         UccAggregator.run(length);
//
//        // Compare all files and write common lines
//        FileComparator.compare(length);
//    }
//}


package at.univie.federated;

import at.univie.enumhyp.GenerateEnumerate;
import at.univie.federated.aggregator.UccAggregator;
import at.univie.federated.client.Client;
import at.univie.federated.compared.FileComparator;
import at.univie.federated.server.Server;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {

    private static String formatDuration(long nanos) {
        long ms = nanos / 1_000_000L;
        long s = ms / 1000L;
        long m = s / 60L;
        long h = m / 60L;
        return String.format("%02dh:%02dm:%02ds (%d ms)", h, (m % 60), (s % 60), ms);
    }

    private static void printTimer(String label, long startNs, long endNs) {
        System.out.println("[TIMER] " + label + ": " + formatDuration(endNs - startNs));
    }

    public static void main(String[] args) throws Exception {

        // ---- Toggle stages here ----
        boolean runEnumhyp = false;
        boolean runServerClient = true;
        boolean runFileCompare = false;
        boolean runAggregator = true;

        String data = "fdReduced";
        int length = 10;

        boolean genGraph = true;
        boolean genEnuTxt = true;

        final long tProgramStart = System.nanoTime();

        // 1) Enumhyp
        if (runEnumhyp) {
            final long tStart = System.nanoTime();

            GenerateEnumerate runGenerateEnumerate = new GenerateEnumerate();
            runGenerateEnumerate.generateEnumerate(data, length, genGraph, genEnuTxt);

            final long tEnd = System.nanoTime();
            printTimer("Enumhyp (generate + enumerate)", tStart, tEnd);
        }

        // 2) Server + Clients
        if (runServerClient) {
            final long tStart = System.nanoTime();

            Server server = new Server();
            Thread serverThread = new Thread(() -> {
                try {
                    server.activateServer();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
            serverThread.setDaemon(true);
            serverThread.start();

            Thread.sleep(1000);

            List<Thread> clientThreads = new ArrayList<>();
            for (int index = 0; index < length; index++) {
                Client client = new Client();
                final int currentPart = index + 1;

                Thread clientThread = new Thread(() -> {
                    client.activateClient(
                            data + "/" + data + "_p" + currentPart + ".txt",
                            "part_" + currentPart
                    );
                });

                clientThreads.add(clientThread);
                clientThread.start();
            }

            for (Thread thread : clientThreads) {
                thread.join();
            }

            final long tEnd = System.nanoTime();
            printTimer("Server + all clients", tStart, tEnd);

            System.out.println("All clients completed. Server continues running...");
        }

        // 3) FileComparator
        if (runFileCompare) {
            final long tStart = System.nanoTime();

            FileComparator.compare(length);

            final long tEnd = System.nanoTime();
            printTimer("FileComparator.compare(" + length + ")", tStart, tEnd);
        }

        // 4) UccAggregator
        if (runAggregator) {
            final long tStart = System.nanoTime();

            UccAggregator.run(length);

            final long tEnd = System.nanoTime();
            printTimer("UccAggregator.run(" + length + ")", tStart, tEnd);
        }

        final long tProgramEnd = System.nanoTime();
        printTimer("TOTAL program", tProgramStart, tProgramEnd);
    }
}
