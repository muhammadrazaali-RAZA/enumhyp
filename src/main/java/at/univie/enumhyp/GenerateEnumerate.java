package at.univie.enumhyp;

import java.io.IOException;
import java.nio.file.Paths;

public final class GenerateEnumerate {

    public void generateEnumerate(String csvFile, int length, boolean graph, boolean txt) throws IOException, InterruptedException {
        EnumhypRunner runner = EnumhypRunner.fromDefaultPath();
        
        // Use csvFile as directory name (e.g., "fdReduced" for data/fdReduced/fdReduced_p1.csv)
        String dir = csvFile;

        // 1. Generate UCC graphs
        if (graph) {
            for (int index = 0; index < length; index++) {
                int currentPart = index + 1;
                runner.generate(dir, Paths.get(csvFile + "_p" + currentPart + ".csv"));
            }
        }

        // 2. Enumerate with output file
        if (txt) {
            for (int index = 0; index < length; index++) {
                int currentPart = index + 1;
                runner.enumerate(
                        dir,
                        Paths.get(csvFile + "_p" + currentPart + "_ucc.graph"),
                        Paths.get(csvFile + "_p" + currentPart + ".txt")
                );
            }
        }
    }


    public void generateEnumerateSingle(String csvFile, boolean graph, boolean txt) throws IOException, InterruptedException {
        // Use centralized path for runner
        EnumhypRunner runner = EnumhypRunner.fromCentralizedPath();

        // Use csvFile as directory name (e.g., "fdReduced" for centralized/data/fdReduced/fdReduced.csv)
        String dir = csvFile;

        // 1. Generate UCC graphs
        if (graph) {
            runner.generate(dir, Paths.get(csvFile + ".csv"));
        }

        // 2. Enumerate with output file
        if (txt) {
            runner.enumerate(
                        dir,
                        Paths.get(csvFile + "_ucc.graph"),
                        Paths.get(csvFile + ".txt")
            );
        }
    }
}