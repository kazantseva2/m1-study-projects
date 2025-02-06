package org.vmk.dep508.library;

import java.sql.*;
import java.util.ArrayList;
import java.util.List;

public class LibraryImpl implements Library {
    private String jdbcUrl;
    private String user;
    private String password;

    public LibraryImpl(String jdbcUrl, String user, String password) {
        this.jdbcUrl = jdbcUrl;
        this.user = user;
        this.password = password;
    }

    private Connection getConnection() throws SQLException {
        return DriverManager.getConnection(jdbcUrl, user, password);
    }

    @Override
    public void addNewBook(Book book) {
        String addSql = "insert into books (book_id, book_title, student_id) values(?, ?, -1)";
        try(PreparedStatement stmp = getConnection().prepareStatement(addSql)){
            stmp.setInt(1, book.getId());
            stmp.setString(2, book.getTitle());
            stmp.execute();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void addAbonent(Student student) {
        String addSql = "insert into abonents (student_id, student_name) values(?, ?)";
        try(PreparedStatement stmp = getConnection().prepareStatement(addSql)){
            stmp.setInt(1, student.getId());
            stmp.setString(2, student.getName());
            stmp.execute();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }


    @Override
    public void borrowBook(Book book, Student student) {
        String sql = "select * from books where student_id = ?";
        try (PreparedStatement stmt = getConnection().prepareStatement(sql)) {
            stmt.setInt(1, student.getId());
            try (ResultSet rs = stmt.executeQuery()) {
                if (!rs.next()) {
                    addAbonent(student);
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }

        String updateSql = "update books set student_id = ? where book_id = ?";
        try(PreparedStatement stmp = getConnection().prepareStatement(updateSql)){
            stmp.setInt(1, student.getId());
            stmp.setInt(2, book.getId());
            stmp.execute();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void returnBook(Book book, Student student) {
        String updateSql = "update books set student_id = -1 where book_id = ?";
        try(PreparedStatement stmp = getConnection().prepareStatement(updateSql)){
            stmp.setInt(1, book.getId());
            stmp.execute();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public List<Book> findAvailableBooks() {
        List<Book> result = new ArrayList<>();
        String sql = "select book_id, book_title from books where student_id = -1";
        try (Statement stmt = getConnection().createStatement();
             ResultSet rs = stmt.executeQuery(sql);) {
            while (rs.next()) {
                Integer bookId = rs.getInt("book_id");
                String bookName = rs.getString("book_title");
                result.add(new Book(bookId, bookName));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return result;
    }

    @Override
    public List<Student> getAllStudents() {
        List<Student> result = new ArrayList<>();
        String sql = "select student_id, student_name from abonents";
        try (Statement stmt = getConnection().createStatement();
             ResultSet rs = stmt.executeQuery(sql);) {
            while (rs.next()) {
                Integer studentId = rs.getInt("student_id");
                String studentName = rs.getString("student_name");
                result.add(new Student(studentId, studentName));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return result;
    }
}
