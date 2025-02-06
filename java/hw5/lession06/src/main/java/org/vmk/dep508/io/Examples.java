package org.vmk.dep508.io;

import org.apache.log4j.Logger;

import java.io.*;


public class Examples {

    private static Logger log = Logger.getLogger(Examples.class);

    public static void main(String[] args) throws IOException {

//        try(InputStream in = new FileInputStream("pom.xml")) {
//            int b;
//            while ((b = in.read()) > 0) {
//                System.out.print((char)b);
//            }
//        }

//        try(InputStream is = new FileInputStream("pom.xml")) {
//            byte[] bytes = new byte[512];
//            while(is.read(bytes) >= 0) {
//                System.out.println(new String(bytes));
//            }
//        } catch (FileNotFoundException e) {
//            e.printStackTrace();
//        }

//        try(FileReader r = new FileReader("pom.xml")) {
//            int b;
//            while ((b = r.read()) >= 0) {
//                System.out.print((char) b);
//            }
//        }


//        try(BufferedReader reader = new BufferedReader(
//                new InputStreamReader(
//                        new FileInputStream("pom.xml")), 1024*1024)) {
//            String line = null;
//            while((line = reader.readLine()) != null) {
//                System.out.println(line);
//            }
//        }

//        try (
//                OutputStream os = new FileOutputStream("pom_copy.xml");
//                InputStream is = new FileInputStream("pom.xml");
//        ) {
//                byte[] bytes = new byte[512];
//                int length = 0;
//                while((length = is.read(bytes)) >= 0) {
//                    os.write(bytes, 0, length);
//                }
//        }


//        StringBuffer target = new StringBuffer("");
//        try (Reader r = new FileReader("pom.xml")) {
//            char[] chars = new char[512];
//            int length = 0;
//            while ( (length = r.read(chars)) >= 0) {
//                target.append(chars, 0, length);
//            }
//        }
//        System.out.println(target);
//
//        try (Writer w = new FileWriter("pom_copy.xml")) {
//            w.append(target);
//        }

//        List<String> lines = Files.readAllLines(FileSystems.getDefault().getPath("pom.xml"));
//        lines.stream().forEach(System.out::println);

//        Scanner s = new Scanner(System.in);
//        s.useDelimiter("\\s+");
//        int i = s.nextInt();
//        double d = s.nextDouble();
//        String str = s.nextLine();
//        System.out.println("" + i + " / " + d + " / " + str);

//        URL url = new URL("http://google.com");
//        URLConnection connection = url.openConnection();
//        try (InputStream yandex = url.openStream()) {
//            int b;
//            while ((b = yandex.read()) >= 0) {
//                System.out.print((char) b);
//            }
//        }



        //string to input stream
//        String test = "One of the greatest strength of the Jackson library is the highly customizable serialization and deserialization process.";
//
//        InputStream in1 = new ByteArrayInputStream(test.getBytes(StandardCharsets.UTF_8));
//
//
//        //read string from input stream using reader
//        BufferedReader reader = new BufferedReader(new InputStreamReader(in1));
//        String redFromStream = reader.readLine();
//        log.info(redFromStream);
    }
}
